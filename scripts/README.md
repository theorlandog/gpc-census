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
