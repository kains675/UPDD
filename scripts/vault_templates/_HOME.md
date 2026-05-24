---
type: moc
title: UPDD Vault — Home
tags:
  - moc
  - home
---

# UPDD Vault — Home

This is the **first thing Claude Code reads** at session start (via SessionStart hook in the UPDD repo). It also serves as the human entrypoint.

## Agent write contract

Claude Code may **write** to:
- `10_Sessions/` — append session logs (one file per session, `YYYY-MM-DD_topic.md`)
- `20_Decisions/` — propose new ADRs (status: `proposed`); user promotes to `accepted`

Claude Code is **read-only** for `30_Lab/`, `40_Knowledge/`, `50_Papers/`, `60_Refs/` unless the user explicitly requests a write.

## Quick links

- Code knowledge: [[MOC_Code]]
- Lab notebook: [[MOC_Lab]]
- Papers & critiques: [[MOC_Papers]]
- Decisions (ADRs): [[MOC_Decisions]]
- UPDD repo: `~/repo` (public, code only)

## Recent sessions

```dataview
TABLE date AS "Date", topics AS "Topics", next AS "Next"
FROM "10_Sessions"
WHERE type = "session"
SORT date DESC
LIMIT 10
```

## Cycles in progress

```dataview
TABLE system, ncaa, cycle, ddg_kcal_mol AS "ΔΔG", sigma_btwn AS "σ_btwn", ci AS "CI", status
FROM "30_Lab"
WHERE type = "experiment" AND (status = "running" OR status = "planned")
SORT cycle DESC, date DESC
```

## Open reviewer critiques

```dataview
TABLE paper, severity, response_due AS "Due", file.link AS "Note"
FROM "50_Papers"
WHERE type = "critique" AND (status = "open" OR status = "in-progress")
SORT severity ASC, response_due ASC
```

## ADRs pending acceptance

```dataview
TABLE date, file.link AS "ADR"
FROM "20_Decisions"
WHERE status = "proposed"
SORT date DESC
```

## Conventions

- File names: `YYYY-MM-DD_kebab-topic.md` for sessions, `ADR-NNNN_kebab-name.md` for decisions
- All notes carry YAML frontmatter with `type:` and `tags:`
- Use `[[wikilinks]]` not Markdown links inside the vault
- `[[#headers]]` link to specific sections
- Numbers in frontmatter unquoted (Dataview parses them as numeric)
