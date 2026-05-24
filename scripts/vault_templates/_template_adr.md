---
adr_id: ADR-XXXX
date: <% tp.date.now("YYYY-MM-DD") %>
status: proposed     # proposed | accepted | superseded | deprecated
supersedes: ""
superseded_by: ""
commit: ""           # git SHA when accepted
schema: ""           # e.g. branched_ddg/0.3 if relevant
tags:
  - decision
---

# ADR-XXXX — <short title>

## Context
What problem prompted this decision. The constraint, the surprise, the bug, the reviewer critique. State the situation, not the answer.

## Decision
What we decided. Concrete. The actual change.

## Consequences
- **Positive**: what improves
- **Negative**: what cost we accept
- **Neutral**: side effects worth noting

## Evidence
- Commit: `<sha>`
- Tests: `tests/test_X.py::test_Y`
- Empirical: numbers, before/after (e.g. σ_btwn 10.78 → 3.72)

## Links
- Session: [[YYYY-MM-DD_topic]]
- Related ADRs: [[ADR-YYYY_name]]
- Code knowledge: [[40_Knowledge/...]]
