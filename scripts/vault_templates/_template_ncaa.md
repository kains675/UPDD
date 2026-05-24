---
type: ncaa
ncaa_id: ""           # 3-letter code, e.g. MTR
ncaa_class: ""        # MTR-class | N-methyl | aromatic-extension | ...
parent_aa: ""         # canonical AA it derives from
parameterization: ""  # GAFF2 | OPLS | hybrid
charge_total: 0       # integer formal charge
charge_engine: ""     # RESP | AM1-BCC | direct
registry_version: ""  # e.g. ncaa_registry/0.7
amber14_patch: ""     # ON | OFF | n/a, q_N value
status: ""            # verified | provisional | deferred
tags:
  - ncaa
---

# <ncaa_id> — <full name>

## Identity
- SMILES: ``
- 3-letter / 1-letter: 
- Parent: 
- Registry entry: `utils/ncaa_registry.py:LINE`

## Parameterization
- Force field: 
- Charge derivation: 
- q_N (backbone N partial charge):  e (AMBER14 patch state)
- Charge consistency audit: `tests/test_charge_consistency.py::test_<ncaa>`

## Structure
```smiles
<SMILES here — rendered by Chem plugin>
```

## Verified systems
- [[<system>_<ncaa>]]

## Known issues
- 

## Links
- ADR: [[ADR-XXXX]]
- Related ncAAs: [[ ]]
