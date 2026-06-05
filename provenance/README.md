# Provenance Log

`runs.jsonl` is an **append-only log** of every run that produced an artifact in
this repo. It is committed to git and is the durable source of truth for figure
lineage.

## Format

Each line is a JSON object with:

| Field | Description |
|---|---|
| `run_id` | UUID for this run |
| `kind` | `figure`, `exploration`, `generic`, or `harvest` |
| `figure` | Figure atom id (for figure runs) |
| `command` | The exact command that was run |
| `operator` | Who ran it |
| `timestamp` | ISO 8601 UTC |
| `outputs` | List of output paths + sha256 hashes |
| `inputs` | Input data paths |
| `repo_sha` | Git commit SHA of this repo at run time |
| `env_hash` | Hash of the Python environment |

## Inspecting

```bash
# Summary view
lab-analysis status --provenance

# Raw log
cat provenance/runs.jsonl | python -m json.tool
```
