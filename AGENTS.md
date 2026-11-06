# Project wiki

`wiki/` is the sole home for project documents: plans, specifications, audits,
design notes, marketing and guides. Keep repository entry points, agent
instructions, licenses, changelogs and files referenced by code, builds, CI,
scripts or other kept files in place. Preserve original documents in `wiki/raw/`;
use topic pages for maintained knowledge and `wiki/log/` for dated activity.

The qmd collection is `deepn-wiki`. The orca registers it against the shared
checkout's `wiki/` after the bootstrap merges; never register a worktree.

After registration, search before reading individual pages:

```bash
qmd search -c deepn-wiki "<literal>" --json -n 5
qmd get "qmd://deepn-wiki/<path>:<from>:<count>"
```

Check status and dates before relying on a page. Use the Check → File → Ingest →
Record → Log loop:

1. **Check:** consult the wiki before starting tasks (`wiki-ask`).
2. **File:** put new documents in the wiki and index originals (`wiki-file`).
3. **Ingest:** integrate sources into topic pages; refresh the registered index
   with `qmd update && qmd embed` after wiki changes, from the shared checkout.
4. **Record:** capture durable decisions and contracts in topic pages.
5. **Log:** record significant work in `wiki/log/YYYY-MM-DD.md` (`wiki-log`).
