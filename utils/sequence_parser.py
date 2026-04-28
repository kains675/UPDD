#!/usr/bin/env python
"""utils/sequence_parser.py — UPDD Sequence Input Parser.

Schema: ``sequence_parser/0.1`` (PR-NEW-A Day 1, 2026-04-28)

Capability Level 1 user-facing input layer (Stage 1 of the v0.7 pipeline,
per `v07_user_verdict/v0.7 Section 4 finalized.md` §4.1.2). Parses user
sequence input in multiple formats and produces a JSON-serializable
canonical representation that downstream stages (PR-NEW-B Branched ΔΔG
Engine via `utils/branched_ddg.py`, PR-NEW-D User CLI) consume.

Supported input formats:
    1. Canonical 1-letter (no ncAA): ``"ICVVQDWGHHRCT"``
    2. Parenthesis-tagged ncAA inline: ``"SLL(MTR)GEYLL"``
       Nested parens for ncAA literal labels: ``"SLL(Trp(1Me))GE"``
    3. Dict-style position map: ``{1: "S", 2: "L", 3: "L", 4: "MTR"}``
       (int or str keys, single-letter values for canonical residues,
       3-letter / synonym values for ncAAs)
    4. PDB-style residue numbering: ``"1:S 2:L 3:L 4:MTR"``

Public API:
    - :class:`SequenceParseResult`  — schema-typed dataclass
    - :class:`SequenceFormatError`  — unrecognized format
    - :class:`NcAAValidationError`  — ncAA code not in registry / unresolved
    - :func:`parse`                  — main entry
    - :func:`detect_format`          — sniff format
    - :func:`resolve_ncaa_synonym`   — synonym → canonical 3-letter code
    - :func:`validate_sequence`      — advisory warnings list

CLI::

    python -m utils.sequence_parser --input "SLL(MTR)GEYLL"
    python -m utils.sequence_parser --input "SLL(1MeW)GE" --resolve-synonyms
    python -m utils.sequence_parser --input '{"1":"S","4":"MTR"}' --format dict

Attribution: synonym lookup table draws on canonical 3-letter codes from
:mod:`utils.ncaa_registry` (UPDD ncAA SSOT, v2). Where the registry already
exposes a synonym (alias tuple), that path is preferred; the inline
:data:`NCAA_SYNONYMS` table acts as a Day-1 supplement for human-readable
spellings (``1MeW``, ``pSer``, ``Trp(1Me)``) commonly seen in ncAA
literature (Khoury 2014, Capece 2012) and user freeform inputs.
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass, field, asdict
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional, Tuple, Union

# ---------------------------------------------------------------------------
# Constants / Schema
# ---------------------------------------------------------------------------

SCHEMA_VERSION = "sequence_parser/0.1"
TOOL_NAME = "utils.sequence_parser"

# Canonical 1-letter → 3-letter map for the 20 standard amino acids.
# Used to (a) validate single-letter tokens, (b) compute formal-charge
# estimate, (c) round-trip dict-style input.
ONE_TO_THREE: Dict[str, str] = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
    "E": "GLU", "Q": "GLN", "G": "GLY", "H": "HIS", "I": "ILE",
    "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
    "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
}
THREE_TO_ONE: Dict[str, str] = {v: k for k, v in ONE_TO_THREE.items()}

# Standard residue formal-charge map (pH 7, default protonation states).
# HIS treated as neutral (HID/HIE tautomer ambiguity → 0). Used by
# :func:`_charge_estimate` for advisory total-charge reporting.
STANDARD_FORMAL_CHARGE: Dict[str, int] = {
    "ARG": +1, "LYS": +1, "HIS": 0,
    "ASP": -1, "GLU": -1,
    "ALA": 0, "ASN": 0, "CYS": 0, "GLN": 0, "GLY": 0, "ILE": 0,
    "LEU": 0, "MET": 0, "PHE": 0, "PRO": 0, "SER": 0, "THR": 0,
    "TRP": 0, "TYR": 0, "VAL": 0,
}

# Day-1 synonym supplement (in addition to ncaa_registry alias tuples).
# Keys are user-facing spellings; values are canonical 3-letter codes from
# the UPDD ncAA SSOT. Lower-cased keys are matched case-insensitively at
# runtime; mixed-case keys preserved here for human readability.
NCAA_SYNONYMS: Dict[str, str] = {
    # 1-Methyl-Tryptophan → MTR
    "1MeW": "MTR", "MeW": "MTR", "1mw": "MTR", "1-MeW": "MTR",
    "Trp(1Me)": "MTR", "1-methyl-Trp": "MTR", "1-MeTrp": "MTR",
    "1-Me-Trp": "MTR", "1-methyltryptophan": "MTR",
    # N-methyl-Alanine → NMA
    "N-MeA": "NMA", "Nme-Ala": "NMA", "MeAla": "NMA", "N-Me-Ala": "NMA",
    # N-methyl-Leucine → NML (registry alias MLE also resolves)
    "N-MeL": "NML", "Nme-Leu": "NML", "MeLeu": "NML", "MLE": "NML",
    # N-methyl-Valine → NMV
    "N-MeV": "NMV", "MeVal": "NMV", "MVA": "NMV",
    # N-methyl-Phenylalanine → MEA
    "N-MeF": "MEA", "MePhe": "MEA", "NMF": "MEA",
    # N-methyl-Glycine (sarcosine) → SAR
    "N-MeG": "SAR", "Sarc": "SAR",
    # D-Alanine → DAL
    "D-Ala": "DAL", "dAla": "DAL",
    # D-Leucine → DLE
    "D-Leu": "DLE", "dLeu": "DLE",
    # D-Tryptophan → DTR
    "D-Trp": "DTR", "dTrp": "DTR",
    # D-Phenylalanine → DPN
    "D-Phe": "DPN", "dPhe": "DPN",
    # D-Tyrosine → DTY
    "D-Tyr": "DTY", "dTyr": "DTY",
    # D-Valine → DVA
    "D-Val": "DVA", "dVal": "DVA",
    # D-Arginine → DAR
    "D-Arg": "DAR", "dArg": "DAR",
    # D-Proline → DPR
    "D-Pro": "DPR", "dPro": "DPR",
    # Hydroxyproline → HYP
    "Hyp": "HYP", "hydroxyproline": "HYP", "trans-Hyp": "HYP",
    # Phospho-Serine → SEP
    "pSer": "SEP", "P-Ser": "SEP", "phospho-Ser": "SEP", "pS": "SEP",
    # Phospho-Threonine → TPO
    "pThr": "TPO", "P-Thr": "TPO", "phospho-Thr": "TPO", "pT": "TPO",
    # Phospho-Tyrosine → PTR
    "pTyr": "PTR", "P-Tyr": "PTR", "phospho-Tyr": "PTR", "pY": "PTR",
    # Norleucine → NLE
    "Nle": "NLE", "norleucine": "NLE",
    # Aminoisobutyric acid → AIB
    "Aib": "AIB", "aminoisobutyric": "AIB",
    # Ornithine → ORN
    "Orn": "ORN", "ornithine": "ORN",
    # Cyclohexylalanine → CHA
    "Cha": "CHA", "Chx-Ala": "CHA",
}


# ---------------------------------------------------------------------------
# Public exception classes
# ---------------------------------------------------------------------------

class SequenceParserError(Exception):
    """Base class for sequence-parser errors."""


class SequenceFormatError(SequenceParserError):
    """Raised when input format cannot be detected or is malformed
    (unbalanced parens, empty token, invalid PDB-style numbering, etc.).
    """


class NcAAValidationError(SequenceParserError):
    """Raised when a token resembling an ncAA code is not resolvable
    against the registry (or, with ``resolve_synonyms=True``, against
    both the registry and :data:`NCAA_SYNONYMS`).

    Carries a ``suggestions`` attribute listing the closest known codes
    (Day-1 implementation: simple shared-prefix scan; full Levenshtein
    fuzzy match deferred to Day-2 per dispatch escalation policy).
    """

    def __init__(self, message: str, *, token: str = "",
                 suggestions: Optional[List[str]] = None):
        super().__init__(message)
        self.token = token
        self.suggestions = suggestions or []


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------

@dataclass
class NcAAPosition:
    """One ncAA occurrence in a parsed sequence."""

    position: int                 # 1-indexed residue position
    code: str                     # canonical 3-letter ncAA code (registry)
    synonym_resolved: bool        # True if the original token was a synonym
    original_token: str           # raw token as the user typed it

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class SequenceParseResult:
    """Schema-typed sequence-parser result.

    Returned by :func:`parse`. Use :meth:`to_dict` for a JSON-serializable
    payload conforming to schema ``sequence_parser/0.1``.

    Fields:
        schema, tool, generated_at: provenance.
        raw_input: the user-supplied input (string form; dict inputs are
            re-serialised to JSON for lossless capture).
        format_detected: ``"canonical" | "parenthesis" | "dict" |
            "pdb_style"``.
        canonical_seq: 1-letter sequence with ncAA positions denoted by
            ``"X"`` (unknown/non-standard placeholder, used widely in
            sequence databases). Length always == ``length``.
        residue_three_letter: full 3-letter residue list (e.g.
            ``["SER","LEU","LEU","MTR","GLY",...]``); preserves ncAA
            identity that ``canonical_seq`` cannot.
        ncaa_positions: list of :class:`NcAAPosition` records.
        length: residue count.
        validation: ``{"status": "PASS"|"WARN", "warnings": [...]}``.
        estimated_total_charge: sum of standard formal charges + ncAA
            formal charges (advisory; HIS counted as neutral).
    """

    schema: str
    tool: str
    generated_at: str
    raw_input: str
    format_detected: str
    canonical_seq: str
    residue_three_letter: List[str]
    ncaa_positions: List[NcAAPosition]
    length: int
    validation: Dict[str, Any]
    estimated_total_charge: int

    def to_dict(self) -> Dict[str, Any]:
        return {
            "schema": self.schema,
            "tool": self.tool,
            "generated_at": self.generated_at,
            "raw_input": self.raw_input,
            "format_detected": self.format_detected,
            "canonical_seq": self.canonical_seq,
            "residue_three_letter": list(self.residue_three_letter),
            "ncaa_positions": [p.to_dict() for p in self.ncaa_positions],
            "length": self.length,
            "validation": dict(self.validation),
            "estimated_total_charge": self.estimated_total_charge,
        }


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

# Cache for the registry index. Lazy-loaded so this module imports cleanly
# even on systems where the optional RDKit / SMILES validation in
# ``ncaa_registry`` errors (we only need code/alias lookup).
_REGISTRY_CACHE: Optional[Dict[str, Any]] = None


def _load_ncaa_registry() -> Dict[str, Any]:
    """Lazy-load the UPDD ncAA SSOT and build a flat lookup index.

    Returns a dict with two keys:
        ``codes``: ``{CODE_UPPER: NCAADef}`` for canonical 3-letter codes.
        ``aliases``: ``{alias_lower: NCAADef}`` covering label/code/PDB/XML
            resnames + every entry in ``NCAADef.aliases``.

    On import failure (e.g. broken environment, missing RDKit forcing
    ``ncaa_registry`` to raise during validation), returns a minimal
    fallback index covering the codes named by the dispatch escalation
    note (MTR, NMA, NML/MLE, DLE, DTR, HYP, SEP, TPO, PTR + a few common
    extensions) so this module stays usable in CI / test sandboxes.
    """
    global _REGISTRY_CACHE
    if _REGISTRY_CACHE is not None:
        return _REGISTRY_CACHE

    out: Dict[str, Any] = {"codes": {}, "aliases": {}}

    try:
        # Defer the import: keep this module light if ncaa_registry has
        # an environment-dependent failure mode.
        from utils import ncaa_registry as _reg  # type: ignore
        for spec in _reg.NCAA_REGISTRY_DATA:
            code_up = spec.code.upper()
            out["codes"][code_up] = spec
            # Alias index
            out["aliases"][spec.code.lower()] = spec
            out["aliases"][spec.label.lower()] = spec
            out["aliases"][spec.pdb_resname.lower()] = spec
            out["aliases"][spec.xml_resname.lower()] = spec
            for a in spec.aliases:
                out["aliases"][a.lower()] = spec
    except Exception:
        # Fallback minimal registry — escalation path per dispatch.
        # These are the codes the dispatch flags as "Day-1 known good".
        _fallback_codes = {
            "MTR": 0, "NMA": 0, "NML": 0, "MLE": 0,
            "DLE": 0, "DTR": 0, "HYP": 0,
            "SEP": -2, "TPO": -2, "PTR": -2,
            "DAL": 0, "DPN": 0, "DAR": +1, "DTY": 0, "DPR": 0, "DVA": 0,
            "NMV": 0, "MEA": 0, "NMK": +1, "NMR": +1, "NMQ": 0, "SAR": 0,
            "AIB": 0, "ORN": +1, "DAB": +1, "HAR": +1, "NLE": 0,
            "TIC": 0, "CHA": 0,
        }

        class _FallbackSpec:
            def __init__(self, code: str, charge: int):
                self.code = code
                self.label = code
                self.pdb_resname = code
                self.xml_resname = code
                self.formal_charge = charge
                self.aliases: Tuple[str, ...] = ()

        for code, ch in _fallback_codes.items():
            spec = _FallbackSpec(code, ch)
            out["codes"][code] = spec
            out["aliases"][code.lower()] = spec

    _REGISTRY_CACHE = out
    return out


def resolve_ncaa_synonym(synonym: str) -> Optional[str]:
    """Return the canonical 3-letter ncAA code for ``synonym``, or ``None``.

    Lookup order:
        1. Exact (case-insensitive) match against :data:`NCAA_SYNONYMS`.
        2. Registry alias index (covers code, label, PDB resname, XML
           resname, and every entry in ``NCAADef.aliases``).

    Returns the canonical code (upper-case) on hit, ``None`` otherwise.
    """
    if not synonym:
        return None
    # Try the inline synonym table first (case-insensitive on keys).
    syn_lc_map = {k.lower(): v for k, v in NCAA_SYNONYMS.items()}
    hit = syn_lc_map.get(synonym.lower())
    if hit is not None:
        return hit.upper()
    # Fall back to the registry alias index.
    reg = _load_ncaa_registry()
    spec = reg["aliases"].get(synonym.lower())
    if spec is not None:
        return spec.code.upper()
    return None


def _suggest_ncaa(token: str, *, max_n: int = 5) -> List[str]:
    """Return up to ``max_n`` registry codes that share a leading
    character or substring with ``token`` (Day-1 best-effort hint;
    fuzzy Levenshtein match deferred per dispatch escalation note).
    """
    if not token:
        return []
    reg = _load_ncaa_registry()
    codes = sorted(reg["codes"].keys())
    tlc = token.lower()
    # Prefer codes starting with the same first character.
    prefix = [c for c in codes if c.lower().startswith(tlc[:1])]
    contain = [c for c in codes if tlc in c.lower() and c not in prefix]
    return (prefix + contain)[:max_n]


def detect_format(raw_input: Union[str, Dict[Any, Any]]) -> str:
    """Sniff the input format. Returns one of:

        ``"dict"``         — a Python dict was passed.
        ``"pdb_style"``    — string with ``N:CODE`` tokens (e.g.
                              ``"1:S 2:L 3:MTR"``).
        ``"parenthesis"``  — string with at least one ``(...)`` block.
        ``"canonical"``    — bare 1-letter sequence string.

    No semantic validation here — see :func:`parse` for that. Empty /
    whitespace-only strings raise :class:`SequenceFormatError`.
    """
    if isinstance(raw_input, dict):
        return "dict"
    if not isinstance(raw_input, str):
        raise SequenceFormatError(
            f"unsupported input type: {type(raw_input).__name__}; "
            "expected str or dict"
        )
    s = raw_input.strip()
    if not s:
        raise SequenceFormatError("empty input string")
    # PDB-style: at least one ``N:TOK`` token (digits, colon, alpha+).
    # Accept whitespace OR comma separators.
    import re as _re
    pdb_token_re = _re.compile(r"\d+\s*:\s*[A-Za-z][A-Za-z0-9_-]*")
    if pdb_token_re.search(s):
        return "pdb_style"
    if "(" in s or ")" in s:
        return "parenthesis"
    return "canonical"


def _split_paren_tokens(seq: str) -> List[Tuple[str, bool]]:
    """Split a parenthesis-tagged sequence into ``(token, is_ncaa)`` pairs.

    Handles nested parens (e.g. ``"SLL(Trp(1Me))GE"`` → ``[("S",False),
    ("L",False),("L",False),("Trp(1Me)",True),("G",False),("E",False)]``).
    Whitespace inside parens is preserved (some labels use it). Whitespace
    outside is stripped.
    """
    tokens: List[Tuple[str, bool]] = []
    i = 0
    n = len(seq)
    while i < n:
        ch = seq[i]
        if ch.isspace():
            i += 1
            continue
        if ch == "(":
            # Find the matching close paren, supporting nesting.
            depth = 1
            j = i + 1
            while j < n and depth > 0:
                if seq[j] == "(":
                    depth += 1
                elif seq[j] == ")":
                    depth -= 1
                j += 1
            if depth != 0:
                raise SequenceFormatError(
                    f"unbalanced parentheses in sequence at offset {i}: "
                    f"{seq!r}"
                )
            inner = seq[i + 1: j - 1].strip()
            if not inner:
                raise SequenceFormatError(
                    f"empty parenthesis token at offset {i}: {seq!r}"
                )
            tokens.append((inner, True))
            i = j
            continue
        if ch == ")":
            raise SequenceFormatError(
                f"stray closing parenthesis at offset {i}: {seq!r}"
            )
        # Single canonical 1-letter token.
        tokens.append((ch, False))
        i += 1
    return tokens


def _resolve_ncaa_token(token: str, *, position: int,
                        resolve_synonyms: bool, strict: bool
                        ) -> Tuple[str, bool, str]:
    """Resolve an ncAA token to ``(canonical_code, synonym_resolved,
    original_token)``.

    Lookup order (uses :func:`resolve_ncaa_synonym`):
        1. Exact registry code (preferred — no synonym flag set).
        2. Synonym table / registry alias (synonym_resolved=True).

    On miss, raises :class:`NcAAValidationError` with suggestions.
    """
    reg = _load_ncaa_registry()
    tok_up = token.upper()
    # Path 1: canonical code (e.g. "MTR")
    if tok_up in reg["codes"]:
        return tok_up, False, token
    # Path 2: synonym resolution (controlled by flag)
    if resolve_synonyms:
        canon = resolve_ncaa_synonym(token)
        if canon is not None:
            return canon, True, token
    # Miss → error
    suggestions = _suggest_ncaa(token)
    if suggestions:
        msg = (
            f"unknown ncAA token {token!r} at position {position}; "
            f"closest registry codes: {', '.join(suggestions)}"
        )
    else:
        msg = (
            f"unknown ncAA token {token!r} at position {position}; "
            f"no close matches in registry"
        )
    raise NcAAValidationError(msg, token=token, suggestions=suggestions)


# --- Format-specific parsers ------------------------------------------------

def _parse_canonical(seq: str) -> Tuple[List[str], List[NcAAPosition]]:
    """1-letter canonical sequence → (3-letter list, [] ncAA positions).

    Whitespace-tolerant. Raises on any non-standard letter.
    """
    out: List[str] = []
    s = "".join(seq.split())
    for i, ch in enumerate(s, 1):
        ch_up = ch.upper()
        if ch_up not in ONE_TO_THREE:
            raise SequenceFormatError(
                f"unknown 1-letter residue {ch!r} at position {i}; "
                f"canonical alphabet: {sorted(ONE_TO_THREE)}"
            )
        out.append(ONE_TO_THREE[ch_up])
    return out, []


def _parse_parenthesis(seq: str, *, resolve_synonyms: bool, strict: bool
                       ) -> Tuple[List[str], List[NcAAPosition]]:
    """Parens-tagged sequence → (3-letter list, ncAA positions)."""
    tokens = _split_paren_tokens(seq)
    out: List[str] = []
    ncaa: List[NcAAPosition] = []
    for idx, (tok, is_ncaa) in enumerate(tokens, 1):
        if not is_ncaa:
            ch_up = tok.upper()
            if ch_up not in ONE_TO_THREE:
                raise SequenceFormatError(
                    f"unknown 1-letter residue {tok!r} at position {idx}"
                )
            out.append(ONE_TO_THREE[ch_up])
        else:
            code, synres, orig = _resolve_ncaa_token(
                tok, position=idx,
                resolve_synonyms=resolve_synonyms, strict=strict,
            )
            out.append(code)
            ncaa.append(NcAAPosition(
                position=idx, code=code,
                synonym_resolved=synres, original_token=orig,
            ))
    return out, ncaa


def _parse_dict(d: Dict[Any, Any], *, resolve_synonyms: bool,
                strict: bool) -> Tuple[List[str], List[NcAAPosition]]:
    """Dict ``{position: residue}`` → (3-letter list, ncAA positions).

    Keys may be int or str(int). Values may be 1-letter (canonical) or
    3-letter (ncAA / canonical 3-letter). Position numbering must be
    contiguous starting at 1; gaps raise :class:`SequenceFormatError`.
    """
    if not d:
        raise SequenceFormatError("empty dict input")
    # Normalize keys to int.
    norm: Dict[int, str] = {}
    for k, v in d.items():
        try:
            ki = int(k)
        except (TypeError, ValueError):
            raise SequenceFormatError(
                f"non-integer position key in dict input: {k!r}"
            )
        if not isinstance(v, str) or not v.strip():
            raise SequenceFormatError(
                f"non-string or empty residue at position {ki}: {v!r}"
            )
        norm[ki] = v.strip()
    # Verify contiguous numbering.
    sorted_keys = sorted(norm)
    if sorted_keys[0] != 1:
        raise SequenceFormatError(
            f"dict input must start at position 1; first key = "
            f"{sorted_keys[0]}"
        )
    for expected, actual in enumerate(sorted_keys, 1):
        if expected != actual:
            raise SequenceFormatError(
                f"dict input has gap at position {expected} "
                f"(found {actual} instead)"
            )
    out: List[str] = []
    ncaa: List[NcAAPosition] = []
    reg = _load_ncaa_registry()
    for pos in sorted_keys:
        v = norm[pos]
        if len(v) == 1:
            ch_up = v.upper()
            if ch_up not in ONE_TO_THREE:
                raise SequenceFormatError(
                    f"unknown 1-letter residue {v!r} at position {pos}"
                )
            out.append(ONE_TO_THREE[ch_up])
            continue
        # Multi-char value: try canonical 3-letter first, then ncAA path.
        v_up = v.upper()
        if v_up in THREE_TO_ONE:
            out.append(v_up)
            continue
        if v_up in reg["codes"]:
            out.append(v_up)
            ncaa.append(NcAAPosition(
                position=pos, code=v_up,
                synonym_resolved=False, original_token=v,
            ))
            continue
        # Synonym path.
        code, synres, orig = _resolve_ncaa_token(
            v, position=pos,
            resolve_synonyms=resolve_synonyms, strict=strict,
        )
        out.append(code)
        ncaa.append(NcAAPosition(
            position=pos, code=code,
            synonym_resolved=synres, original_token=orig,
        ))
    return out, ncaa


def _parse_pdb_style(s: str, *, resolve_synonyms: bool, strict: bool
                     ) -> Tuple[List[str], List[NcAAPosition]]:
    """PDB-style ``"N:CODE [N:CODE ...]"`` → (3-letter list, ncAA positions).

    Tokens may be separated by whitespace OR commas. Internally builds a
    dict-input and dispatches to :func:`_parse_dict` for shared validation.
    """
    import re as _re
    s = s.strip()
    # Split on whitespace and/or commas.
    chunks = [c for c in _re.split(r"[\s,]+", s) if c]
    d: Dict[int, str] = {}
    for ch in chunks:
        if ":" not in ch:
            raise SequenceFormatError(
                f"PDB-style token missing ':' separator: {ch!r}"
            )
        pos_s, val = ch.split(":", 1)
        pos_s = pos_s.strip()
        val = val.strip()
        try:
            pos = int(pos_s)
        except ValueError:
            raise SequenceFormatError(
                f"PDB-style position must be an integer: {pos_s!r}"
            )
        if not val:
            raise SequenceFormatError(
                f"PDB-style residue must be non-empty: {ch!r}"
            )
        if pos in d:
            raise SequenceFormatError(
                f"duplicate PDB-style position: {pos}"
            )
        d[pos] = val
    return _parse_dict(d, resolve_synonyms=resolve_synonyms, strict=strict)


# ---------------------------------------------------------------------------
# Charge estimate / validation
# ---------------------------------------------------------------------------

def _charge_estimate(residue_three_letter: List[str]) -> int:
    """Sum of standard formal charges + ncAA registry charges (advisory).

    Standard residues: per :data:`STANDARD_FORMAL_CHARGE` (HIS=0).
    ncAA residues:    per ``NCAADef.formal_charge`` from the registry,
                      else 0 if the code is unknown (best-effort; the
                      caller has already raised on unresolved tokens by
                      the time we get here).
    """
    reg = _load_ncaa_registry()
    total = 0
    for code in residue_three_letter:
        cu = code.upper()
        if cu in STANDARD_FORMAL_CHARGE:
            total += STANDARD_FORMAL_CHARGE[cu]
            continue
        spec = reg["codes"].get(cu)
        if spec is not None:
            total += int(getattr(spec, "formal_charge", 0))
    return total


def validate_sequence(result: SequenceParseResult) -> List[str]:
    """Advisory validation. Returns a list of warning strings (may be
    empty). Does NOT raise — sign that the parse itself succeeded.

    Checks:
        - Empty sequence (length == 0)               → "empty sequence"
        - Length > 100                                → "long sequence"
        - |estimated_total_charge| > 6                → "unusual charge"
        - All-ncAA sequence (no canonical residues)   → advisory only
    """
    warnings: List[str] = []
    if result.length == 0:
        warnings.append("empty sequence (length 0)")
    if result.length > 100:
        warnings.append(
            f"long sequence (length {result.length} > 100); "
            "verify this is intentional"
        )
    if abs(result.estimated_total_charge) > 6:
        warnings.append(
            f"unusual estimated total charge "
            f"{result.estimated_total_charge:+d} "
            "(|q| > 6 is rare for short peptides)"
        )
    if result.length > 0 and len(result.ncaa_positions) == result.length:
        warnings.append(
            "all-ncAA sequence (no canonical residues); confirm intent"
        )
    return warnings


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def parse(raw_input: Union[str, Dict[Any, Any]], *,
          resolve_synonyms: bool = True,
          strict: bool = False,
          format_hint: Optional[str] = None) -> SequenceParseResult:
    """Parse ``raw_input`` into a :class:`SequenceParseResult`.

    Args:
        raw_input: User input (string or dict — see module docstring).
        resolve_synonyms: If True (default), unrecognised ncAA tokens are
            also looked up in :data:`NCAA_SYNONYMS` and the registry alias
            index. If False, only canonical 3-letter codes match.
        strict: Reserved. Currently a no-op pass-through; populated for
            forward compatibility with stricter validation modes (e.g.
            forbidding all-ncAA sequences).
        format_hint: If provided, skip :func:`detect_format` and force the
            given format. One of ``"canonical" | "parenthesis" | "dict" |
            "pdb_style" | "auto"``.

    Returns:
        SequenceParseResult.

    Raises:
        SequenceFormatError: malformed input.
        NcAAValidationError: unresolved ncAA token (may carry suggestions).
    """
    if format_hint and format_hint != "auto":
        fmt = format_hint
        # Validate type-format compatibility for dict.
        if fmt == "dict" and not isinstance(raw_input, dict):
            # Best-effort: try parsing as JSON.
            try:
                parsed = json.loads(raw_input) if isinstance(raw_input, str) else None
            except (json.JSONDecodeError, TypeError):
                parsed = None
            if isinstance(parsed, dict):
                raw_input = parsed
            else:
                raise SequenceFormatError(
                    "format_hint='dict' requires a dict input or JSON-encoded "
                    "object string"
                )
    else:
        fmt = detect_format(raw_input)

    # Capture raw_input as a string for provenance (round-trip stable).
    if isinstance(raw_input, dict):
        raw_str = json.dumps(raw_input, sort_keys=True)
    else:
        raw_str = str(raw_input)

    if fmt == "canonical":
        residues, ncaas = _parse_canonical(str(raw_input))
    elif fmt == "parenthesis":
        residues, ncaas = _parse_parenthesis(
            str(raw_input), resolve_synonyms=resolve_synonyms, strict=strict
        )
    elif fmt == "dict":
        d = raw_input if isinstance(raw_input, dict) else json.loads(str(raw_input))
        residues, ncaas = _parse_dict(
            d, resolve_synonyms=resolve_synonyms, strict=strict
        )
    elif fmt == "pdb_style":
        residues, ncaas = _parse_pdb_style(
            str(raw_input), resolve_synonyms=resolve_synonyms, strict=strict
        )
    else:
        raise SequenceFormatError(
            f"unknown format_hint {fmt!r}; expected one of "
            "{canonical, parenthesis, dict, pdb_style, auto}"
        )

    # Build canonical 1-letter sequence (X for ncAAs).
    canonical_chars: List[str] = []
    ncaa_pos_set = {p.position for p in ncaas}
    for i, code in enumerate(residues, 1):
        if i in ncaa_pos_set or code not in THREE_TO_ONE:
            canonical_chars.append("X")
        else:
            canonical_chars.append(THREE_TO_ONE[code])
    canonical_seq = "".join(canonical_chars)

    total_charge = _charge_estimate(residues)

    result = SequenceParseResult(
        schema=SCHEMA_VERSION,
        tool=TOOL_NAME,
        generated_at=datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        raw_input=raw_str,
        format_detected=fmt,
        canonical_seq=canonical_seq,
        residue_three_letter=residues,
        ncaa_positions=ncaas,
        length=len(residues),
        validation={"status": "PASS", "warnings": []},
        estimated_total_charge=total_charge,
    )

    warnings = validate_sequence(result)
    if warnings:
        result.validation = {"status": "WARN", "warnings": warnings}

    return result


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="python -m utils.sequence_parser",
        description=(
            f"UPDD Sequence Input Parser (PR-NEW-A; schema {SCHEMA_VERSION})."
        ),
    )
    ap.add_argument("--input", required=True,
                    help="Sequence input string or JSON-encoded dict.")
    ap.add_argument("--resolve-synonyms", action="store_true", default=True,
                    help="Resolve ncAA synonyms (default ON).")
    ap.add_argument("--no-resolve-synonyms", dest="resolve_synonyms",
                    action="store_false",
                    help="Disable synonym resolution; require canonical codes.")
    ap.add_argument("--strict", action="store_true",
                    help="Reserved; tightens validation in future revs.")
    ap.add_argument("--format", choices=("auto", "canonical", "parenthesis",
                                         "parens", "dict", "pdb_style", "pdb"),
                    default="auto",
                    help="Force a format. Default 'auto' = detect_format().")
    ap.add_argument("--output", default=None,
                    help="If set, write JSON payload to this path. Stdout "
                         "always carries a one-line summary.")
    return ap


def main(argv: Optional[List[str]] = None) -> int:
    ap = _build_parser()
    args = ap.parse_args(argv)
    fmt_norm = {
        "auto": "auto",
        "canonical": "canonical",
        "parenthesis": "parenthesis", "parens": "parenthesis",
        "dict": "dict",
        "pdb_style": "pdb_style", "pdb": "pdb_style",
    }[args.format]
    try:
        result = parse(
            args.input,
            resolve_synonyms=args.resolve_synonyms,
            strict=args.strict,
            format_hint=fmt_norm,
        )
    except NcAAValidationError as exc:
        print(f"[sequence_parser] NcAAValidationError: {exc}", file=sys.stderr)
        if exc.suggestions:
            print(
                f"[sequence_parser]   suggestions: {', '.join(exc.suggestions)}",
                file=sys.stderr,
            )
        return 4
    except SequenceFormatError as exc:
        print(f"[sequence_parser] SequenceFormatError: {exc}", file=sys.stderr)
        return 3
    except SequenceParserError as exc:
        print(f"[sequence_parser] ERROR: {type(exc).__name__}: {exc}",
              file=sys.stderr)
        return 2

    payload = result.to_dict()
    print(f"[sequence_parser] schema={result.schema}")
    print(f"[sequence_parser] format={result.format_detected}  "
          f"length={result.length}  ncaas={len(result.ncaa_positions)}  "
          f"q={result.estimated_total_charge:+d}")
    print(f"[sequence_parser] canonical_seq={result.canonical_seq}")
    if result.ncaa_positions:
        ncaa_summary = ", ".join(
            f"{p.position}:{p.code}"
            + ("(syn)" if p.synonym_resolved else "")
            for p in result.ncaa_positions
        )
        print(f"[sequence_parser] ncaa_positions=[{ncaa_summary}]")
    print(f"[sequence_parser] validation={result.validation['status']}")
    for w in result.validation.get("warnings", []):
        print(f"[sequence_parser]   warning: {w}")

    if args.output:
        with open(args.output, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)
        print(f"[sequence_parser] wrote: {args.output}")
    else:
        # When no --output, also dump JSON to stdout (after summary lines)
        # so this CLI can be piped (`python -m utils.sequence_parser ... | jq`).
        print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
