#!/usr/bin/env python
"""
D1 QM stack cross-validation: PySCF (qmmm env) vs Psi4 (psi4 env) vs Khoury published.

Loads the most-recent output of:
  - scripts/d1_khoury_resp_reproduction.py  (PySCF stack)
  - scripts/d1_khoury_psi4_resp.py          (Psi4 stack)

and generates a 3-way per-atom comparison + 3 pairwise verdicts (PySCF vs
Khoury, Psi4 vs Khoury, Psi4 vs PySCF).

Reproduction PASS = max |Δq| < 0.005 e per atom in the relevant pairing.

Output: outputs/analysis/d1_qmstack_crossvalidation_{TS}/
  - comparison.md       (3-way report)
  - comparison.json     (machine-readable summary)

Run:
  python scripts/d1_qmstack_crossvalidation.py
"""
from __future__ import annotations

import json
import sys
from datetime import datetime
from pathlib import Path

PROJ = Path(__file__).resolve().parent.parent

KHOURY_OMW_CHARGES = {
    "N":   -0.280943, "H":   0.241023, "CA": -0.124259, "HA":  0.120930,
    "C":    0.575115, "O":  -0.571038, "CB": -0.157536, "HB2": 0.104569,
    "HB3":  0.104569, "CG":  -0.126422, "CD2": 0.083375, "CE2": 0.119586,
    "NE1":  0.077670, "CD1": -0.252655, "HD1": 0.204947, "CM": -0.239062,
    "HM1":  0.112200, "HM2":  0.112200, "HM3":  0.112200, "CZ2": -0.241414,
    "HZ2":  0.146305, "CH2": -0.139177, "HH2":  0.136546, "CZ3": -0.206367,
    "HZ3":  0.143250, "CE3": -0.211658, "HE3":  0.156047,
}

PASS_THRESHOLD = 0.005  # |Δq| per atom


def find_latest(pattern: str, charges_filename: str) -> tuple[Path | None, dict[str, float] | None]:
    candidates = sorted(
        (PROJ / "outputs/analysis").glob(pattern), reverse=True
    )
    for d in candidates:
        f = d / charges_filename
        if f.exists():
            try:
                charges = json.loads(f.read_text())
                charges = {k: float(v) for k, v in charges.items()}
                return d, charges
            except Exception as e:
                print(f"[warn] Failed to load {f}: {e}")
                continue
    return None, None


def pairwise_stats(a: dict[str, float], b: dict[str, float], atom_order: list[str]) -> dict:
    deltas = {}
    sqsum = 0.0
    n = 0
    max_abs = 0.0
    max_atom = None
    for name in atom_order:
        if name not in a or name not in b:
            continue
        d = a[name] - b[name]
        deltas[name] = d
        sqsum += d ** 2
        n += 1
        if abs(d) > max_abs:
            max_abs = abs(d)
            max_atom = name
    rmsd = (sqsum / n) ** 0.5 if n else float("nan")
    verdict = "PASS" if max_abs < PASS_THRESHOLD else "FAIL"
    return {
        "rmsd": rmsd,
        "max_abs": max_abs,
        "max_atom": max_atom,
        "n": n,
        "verdict": verdict,
        "deltas": deltas,
    }


def main() -> int:
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    outdir = PROJ / f"outputs/analysis/d1_qmstack_crossvalidation_{ts}"
    outdir.mkdir(parents=True, exist_ok=True)

    pyscf_dir, pyscf_q = find_latest(
        "d1_khoury_resp_reproduction_*", "omw_pyscf_charges.json"
    )
    psi4_dir, psi4_q = find_latest(
        "d1_psi4_resp_reproduction_*", "omw_psi4_charges.json"
    )

    if pyscf_q is None:
        print(f"[warn] No PySCF charges found (expected omw_pyscf_charges.json in outputs/analysis/d1_khoury_resp_reproduction_*/)")
    if psi4_q is None:
        print(f"[warn] No Psi4 charges found (expected omw_psi4_charges.json in outputs/analysis/d1_psi4_resp_reproduction_*/)")

    if pyscf_q is None and psi4_q is None:
        sys.exit("[error] Neither reproduction output found. Run scripts/d1_khoury_resp_reproduction.py and/or scripts/d1_khoury_psi4_resp.py first.")

    atom_order = list(KHOURY_OMW_CHARGES.keys())

    stats = {}
    if pyscf_q is not None:
        stats["pyscf_vs_khoury"] = pairwise_stats(pyscf_q, KHOURY_OMW_CHARGES, atom_order)
    if psi4_q is not None:
        stats["psi4_vs_khoury"] = pairwise_stats(psi4_q, KHOURY_OMW_CHARGES, atom_order)
    if pyscf_q is not None and psi4_q is not None:
        stats["psi4_vs_pyscf"] = pairwise_stats(psi4_q, pyscf_q, atom_order)

    md = [
        "# D1 QM Stack Cross-Validation Report",
        "",
        f"**Date**: {datetime.now().isoformat()}",
        f"**Goal**: Confirm Khoury 2014 OMW RESP-A2 charges are reproducible across two independent QM stacks (PySCF/Antechamber vs psi4/resp-plugin), establishing stack-agnostic equivalence for v0.9 Layer 3D roadmap.",
        f"**Tolerance**: max |Δq| < {PASS_THRESHOLD} e per atom (per pairing)",
        "",
        "## Sources",
        "",
        f"- **PySCF stack**: `{pyscf_dir.name if pyscf_dir else 'NOT FOUND'}`",
        f"- **Psi4 stack**:  `{psi4_dir.name if psi4_dir else 'NOT FOUND'}`",
        f"- **Khoury reference**: hardcoded from `Reference/Khoury 2014/sb400168u_si_002/ffncaa.in` L778-830",
        "",
        "## Summary verdicts",
        "",
        "| Comparison | n | RMSD (e) | max |Δq| (e) | max-Δ atom | Verdict |",
        "|---|---|---|---|---|---|",
    ]
    for key, label in [
        ("pyscf_vs_khoury", "PySCF vs Khoury"),
        ("psi4_vs_khoury", "Psi4 vs Khoury"),
        ("psi4_vs_pyscf", "Psi4 vs PySCF"),
    ]:
        s = stats.get(key)
        if s is None:
            md.append(f"| {label} | n/a | n/a | n/a | n/a | **MISSING** |")
        else:
            md.append(
                f"| {label} | {s['n']} | {s['rmsd']:.4f} | {s['max_abs']:.4f} | "
                f"{s['max_atom']} | **{s['verdict']}** |"
            )

    md.extend([
        "",
        "## Per-atom 3-way comparison",
        "",
        "| Atom | Khoury (e) | PySCF (e) | Δ(PySCF−Kh) | Psi4 (e) | Δ(Psi4−Kh) | Δ(Psi4−PySCF) |",
        "|---|---|---|---|---|---|---|",
    ])
    for name in atom_order:
        kq = KHOURY_OMW_CHARGES[name]
        pq = pyscf_q.get(name) if pyscf_q else None
        sq = psi4_q.get(name) if psi4_q else None
        d_pq = pq - kq if pq is not None else None
        d_sq = sq - kq if sq is not None else None
        d_inter = (sq - pq) if (pq is not None and sq is not None) else None
        cells = [
            name,
            f"{kq:+.4f}",
            f"{pq:+.4f}" if pq is not None else "—",
            f"{d_pq:+.4f}" if d_pq is not None else "—",
            f"{sq:+.4f}" if sq is not None else "—",
            f"{d_sq:+.4f}" if d_sq is not None else "—",
            f"{d_inter:+.4f}" if d_inter is not None else "—",
        ]
        md.append("| " + " | ".join(cells) + " |")

    md.extend([
        "",
        "## Decision matrix",
        "",
        "| Outcome | Implication for v0.9 Layer 3D |",
        "|---|---|",
        "| **Both PASS vs Khoury** (PySCF & Psi4) | Stack-agnostic reproduction confirmed. D2/D3 may proceed with either stack; recommend psi4 (native RESP plugin, fewer Antechamber bridge steps). |",
        "| **Both FAIL but Psi4 ≈ PySCF** (Δq inter < 0.005) | Systematic stack-agnostic deviation from Khoury — investigate Khoury 2014 SI for non-standard protocol detail (capping convention, ESP grid weighting, multi-conformer fit). |",
        "| **One PASS, one FAIL** | Stack-specific bias — adopt the PASSing stack as Layer 3D reference; document the failing stack's bias source. |",
        "| **Both FAIL and Psi4 ≠ PySCF** | Independent issues per stack — re-examine RESP plugin parameters and atom equivalence groups in both. |",
        "",
        "## Source files",
        "",
        f"- PySCF charges: `{pyscf_dir / 'omw_pyscf_charges.json' if pyscf_dir else 'n/a'}`",
        f"- Psi4 charges:  `{psi4_dir / 'omw_psi4_charges.json' if psi4_dir else 'n/a'}`",
        f"- This report:   `{outdir / 'comparison.md'}`",
        f"- JSON summary:  `{outdir / 'comparison.json'}`",
    ])

    (outdir / "comparison.md").write_text("\n".join(md))

    # JSON summary
    json_summary = {
        "date": datetime.now().isoformat(),
        "pass_threshold_e": PASS_THRESHOLD,
        "sources": {
            "pyscf_dir": str(pyscf_dir) if pyscf_dir else None,
            "psi4_dir": str(psi4_dir) if psi4_dir else None,
        },
        "verdicts": {k: {kk: vv for kk, vv in v.items() if kk != "deltas"}
                     for k, v in stats.items()},
    }
    (outdir / "comparison.json").write_text(json.dumps(json_summary, indent=2))

    print(f"\nReport: {outdir / 'comparison.md'}")
    print(f"Summary JSON: {outdir / 'comparison.json'}")
    for k, s in stats.items():
        print(f"  {k}: max |Δq|={s['max_abs']:.4f} ({s['max_atom']}) → {s['verdict']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
