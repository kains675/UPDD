#!/usr/bin/env bash
# bootstrap_vault.sh — Create ~/UPDD-vault Obsidian vault scaffold.
#
# Idempotent: existing files are never overwritten. Re-run after editing
# scripts/vault_templates/ to add new template files without touching
# anything you've already written in the vault.
#
# Usage:
#   bash scripts/bootstrap_vault.sh                  # default ~/UPDD-vault
#   VAULT_DIR=/path/to/vault bash scripts/bootstrap_vault.sh

set -euo pipefail

VAULT_DIR="${VAULT_DIR:-$HOME/UPDD-vault}"
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TEMPLATE_SRC="$REPO_DIR/scripts/vault_templates"

if [[ ! -d "$TEMPLATE_SRC" ]]; then
    echo "ERROR: template source not found: $TEMPLATE_SRC" >&2
    exit 1
fi

echo "Bootstrapping Obsidian vault at: $VAULT_DIR"
mkdir -p "$VAULT_DIR"

# --- folder skeleton (7 top-level dirs, flat by design) ---
declare -a DIRS=(
    "00_Index"
    "10_Sessions"
    "20_Decisions"
    "30_Lab/cycles"
    "30_Lab/systems"
    "40_Knowledge/ncaa"
    "40_Knowledge/pbc"
    "40_Knowledge/schema"
    "40_Knowledge/statistics"
    "40_Knowledge/pipeline"
    "50_Papers/paper-1-v2"
    "60_Refs/zotero-imports"
    "60_Refs/pdf-notes"
    ".obsidian"
)
for d in "${DIRS[@]}"; do
    mkdir -p "$VAULT_DIR/$d"
done

# --- copy templates (never overwrite) ---
copy_if_absent() {
    local src="$1" dst="$2"
    if [[ -e "$dst" ]]; then
        echo "  skip (exists): ${dst#$VAULT_DIR/}"
    else
        cp "$src" "$dst"
        echo "  wrote:        ${dst#$VAULT_DIR/}"
    fi
}

copy_if_absent "$TEMPLATE_SRC/_HOME.md"               "$VAULT_DIR/00_Index/_HOME.md"
copy_if_absent "$TEMPLATE_SRC/_template_session.md"   "$VAULT_DIR/10_Sessions/_template_session.md"
copy_if_absent "$TEMPLATE_SRC/_template_adr.md"       "$VAULT_DIR/20_Decisions/_template_adr.md"
copy_if_absent "$TEMPLATE_SRC/_template_experiment.md" "$VAULT_DIR/30_Lab/_template_experiment.md"
copy_if_absent "$TEMPLATE_SRC/_template_ncaa.md"      "$VAULT_DIR/40_Knowledge/ncaa/_template_ncaa.md"
copy_if_absent "$TEMPLATE_SRC/_template_critique.md"  "$VAULT_DIR/50_Papers/_template_critique.md"

# --- stub MOC files in 00_Index/ ---
write_if_absent() {
    local path="$1" body="$2"
    if [[ -e "$path" ]]; then
        echo "  skip (exists): ${path#$VAULT_DIR/}"
    else
        printf '%s\n' "$body" > "$path"
        echo "  wrote:        ${path#$VAULT_DIR/}"
    fi
}

write_if_absent "$VAULT_DIR/00_Index/MOC_Code.md" '---
type: moc
title: Code Domain Knowledge
tags: [moc, code]
---

# Code Domain Knowledge

Atomic notes in `40_Knowledge/` organized by area.

## ncAA registry
```dataview
LIST FROM "40_Knowledge/ncaa" WHERE type = "ncaa" SORT file.name ASC
```

## PBC verification layers
```dataview
LIST FROM "40_Knowledge/pbc" SORT file.name ASC
```

## Schemas
```dataview
LIST FROM "40_Knowledge/schema" SORT file.name ASC
```

## Statistics framework
σ_btwn, σ_w, Convergence Index, z_SE
```dataview
LIST FROM "40_Knowledge/statistics" SORT file.name ASC
```

## Pipeline stages
RFdiffusion → ProteinMPNN → AF2 → MD → QM/MM → MM-PBSA
```dataview
LIST FROM "40_Knowledge/pipeline" SORT file.name ASC
```'

write_if_absent "$VAULT_DIR/00_Index/MOC_Lab.md" '---
type: moc
title: Lab Notebook
tags: [moc, lab]
---

# Lab Notebook

## All experiments — by cycle
```dataview
TABLE system, ncaa, status, ddg_kcal_mol AS "ΔΔG", sigma_btwn AS "σ_btwn", ci AS "CI"
FROM "30_Lab"
WHERE type = "experiment"
SORT cycle ASC, date DESC
```

## Systems
```dataview
LIST FROM "30_Lab/systems" SORT file.name ASC
```

## Tier distribution
```dataview
TABLE length(rows) AS "count"
FROM "30_Lab"
WHERE type = "experiment" AND tier
GROUP BY tier
SORT tier ASC
```'

write_if_absent "$VAULT_DIR/00_Index/MOC_Papers.md" '---
type: moc
title: Papers & Reviewer Critiques
tags: [moc, paper]
---

# Papers & Reviewer Critiques

## Open critiques
```dataview
TABLE paper, reviewer, severity, status, response_due AS "Due"
FROM "50_Papers"
WHERE type = "critique"
SORT status ASC, severity ASC
```

## Manuscripts
```dataview
LIST FROM "50_Papers" WHERE file.name = "manuscript" OR file.name = "revision-arc"
```'

write_if_absent "$VAULT_DIR/00_Index/MOC_Decisions.md" '---
type: moc
title: Architecture Decision Records
tags: [moc, decision]
---

# ADR Index

## By status
```dataview
TABLE adr_id, date, schema, commit
FROM "20_Decisions"
WHERE type != "moc"
SORT status ASC, adr_id DESC
GROUP BY status
```

## All ADRs
```dataview
TABLE adr_id, date, status
FROM "20_Decisions"
WHERE adr_id
SORT adr_id DESC
```'

# --- vault .gitignore ---
write_if_absent "$VAULT_DIR/.gitignore" '# Obsidian device-specific state
.obsidian/workspace.json
.obsidian/workspace-mobile.json
.obsidian/cache
.trash/
# Obsidian Git plugin conflict markers
conflict-files-obsidian-git.md
# Local-only notes (never push, even to private remote)
*.local.md
# OS cruft
.DS_Store
Thumbs.db'

# --- .obsidian/community-plugins.json seed (Obsidian reads this on next launch) ---
write_if_absent "$VAULT_DIR/.obsidian/community-plugins.json" '[
  "obsidian-git",
  "dataview",
  "templater-obsidian",
  "obsidian-tasks-plugin",
  "obsidian-kanban",
  "obsidian-chem"
]'

# --- .obsidian/app.json: enable folder-aware features ---
write_if_absent "$VAULT_DIR/.obsidian/app.json" '{
  "alwaysUpdateLinks": true,
  "newLinkFormat": "shortest",
  "useMarkdownLinks": false,
  "attachmentFolderPath": "60_Refs"
}'

# --- git init (only if not already a repo) ---
# We intentionally do NOT auto-commit: the user's vault repo will likely need
# their own signing config / identity. The bootstrap prints the exact commit
# command in the "Next steps" so they run it once they're satisfied.
if [[ ! -d "$VAULT_DIR/.git" ]]; then
    git -C "$VAULT_DIR" init -q -b main
    echo "  git:          initialized empty repo on 'main' (no commit yet)"
else
    echo "  git:          existing repo, skipped init"
fi

cat <<EOF

Vault ready at: $VAULT_DIR

Next steps (manual):
  1. First commit in the vault repo:
        git -C $VAULT_DIR add .
        git -C $VAULT_DIR commit -m "Bootstrap UPDD vault skeleton"
  2. Open in Obsidian: File → Open vault → $VAULT_DIR
  3. Settings → Community plugins → Turn on, then Browse:
     Obsidian Git, Dataview, Templater, Tasks, Kanban, Chem
  4. Templater settings: Template folder location = (any path containing _template_*.md)
  5. Obsidian Git: set Auto commit-and-sync interval = 15 minutes
  6. Create private GitHub repo and push:
        gh repo create UPDD-vault --private --source=$VAULT_DIR --remote=origin --push
  7. Verify Claude Code SessionStart hook (in UPDD repo):
        Start a new \`claude\` session — first message should reference _HOME.md content.
EOF
