# Research pipeline scripts

These scripts are not part of the gpc-census package and are excluded from its
dependency tree on purpose. They require solver libraries pinned as follows:

- census_engine.py: ortools==9.15.6755, pulp
- state_constructor.py: numpy, scipy
- interference8_1.py: ortools==9.15.6755

Install into a scratch environment, for example:
uv venv scratch && uv pip install --python scratch ortools==9.15.6755 pulp numpy scipy

Verdicts in results/data/census were produced by census_engine.py under the
pinned ortools version. See results/data/PROVENANCE.md for the full
certification chain.

## Checkpointing (spot-resilient runs)

The long campaigns (`solve_all.py`, plus `census_engine.py file` and
`state_constructor.py`'s file mode) share `checkpoint.py`, which periodically
snapshots the output file so a spot-instance kill loses at most one interval of
work, and restores it on a fresh instance so a redeploy resumes where the last
one stopped (the campaigns skip already-recorded work once the file is back).

Common flags:

- `--out PATH` persist results (one JSON record per line) and make the run
  resumable. `solve_all.py` always persists (to
  results/data/states_interference.jsonl by default); for the other two, the
  checkpoint machinery only engages once `--out` is given.
- `--checkpoint-interval SECONDS` wall-clock cadence (default 300).
- `--s3-bucket NAME` also mirror each checkpoint to S3 (durable across losing
  the instance's ephemeral disk); requires `boto3` (`uv pip install .[s3]`).
- `--s3-key KEY` object key (default: the output file's name).
- `--s3-profile NAME` AWS profile; omit to use the default credential chain
  (instance role, environment variables, shared config, ...).

On SIGTERM/SIGINT (a spot reclaim sends SIGTERM ~2 min ahead) the scripts flush
a final checkpoint before exiting. Examples:

    uv run scripts/solve_all.py --checkpoint-interval 120 \
      --s3-bucket my-bucket --s3-profile research
    python scripts/census_engine.py file vertices.json \
      --out census.jsonl --s3-bucket my-bucket
