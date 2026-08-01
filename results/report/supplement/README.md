# Supplement for Paper 1

This archive contains the data and verification programs used for the state
preimage and finite-classification results in the accompanying manuscript. It
is intentionally limited to Paper 1. Fiber dynamics, Hamiltonian studies, and
higher-rank exploratory data are not included.

## Contents

- `data/states.jsonl`: 799 exact state records, one per tabulated vertex.
- `data/vertices/`: the nine rational vertex tables.
- `data/constraints.json`: normalized rational H-representations.
- `data/interference_certificates.json`: 12 finite no-root-distinct proofs.
- `data/constraint_source_map.json`: ordered page and row concordance for the
  rank-9 and rank-10 tables in Altunbulak's dissertation.
- `data/source_rows.json`: the source transcription used by the concordance
  checker.
- `scripts/`: separately implemented state, no-design, H-to-V, source-table,
  and cross-link verifiers.
- `SOURCE.json`: source commit and tool versions recorded by the archive
  builder.
- `SHA256SUMS`: a complete manifest for every other file in the archive.

The JSON state records use zero-based orbital indices; formulas in the paper
use one-based indices. Each record contains a system, zero-based vertex index,
canonical integer spectrum and denominator, classification label, and an exact
`closed_form` object. Internal timings, solver tiers, floating discovery
amplitudes, and transport notes are excluded. The exact object contains the
common squared-amplitude denominator, squared-amplitude weights or algebraic
expressions, display expressions for the amplitudes, and the determinant
support.

## Reproduction

Python 3.12.13 and `uv` are recommended. Keep the generated environment outside
the checksummed archive, then create the locked symbolic environment:

```sh
export UV_PROJECT_ENVIRONMENT=/tmp/gpc-census-paper1-venv
uv sync --locked --no-install-project
```

Run the complete cross-link and checksum check:

```sh
uv run --no-sync python scripts/verify_paper1_bundle_standalone.py .
```

Rebuild and verify every one-body density matrix:

```sh
uv run --no-sync python scripts/verify_states_standalone.py data/states.jsonl
```

The finite no-design proofs and the H-to-V reconstruction require only the
Python standard library:

```sh
python3 scripts/verify_interference_certificates_standalone.py \
  data/interference_certificates.json
python3 scripts/verify_vertex_exhaustion_standalone.py \
  --constraints data/constraints.json --vertices data/vertices
```

The dissertation PDF is not redistributed. After downloading the bitstream
identified in `data/constraint_source_map.json`, reproduce the page map with:

```sh
python3 scripts/verify_constraint_source_standalone.py dissertation.pdf \
  --source-rows data/source_rows.json \
  --normalized data/constraints.json \
  --expected-map data/constraint_source_map.json
```

The expected map records the PDF hash and the Poppler `pdftotext` version used
to produce it. The source checker confirms all 60 four-fermion, nine-orbital
rows and all 379 rank-10 rows in order. Appendix A prints only 51 of the stated
52 three-fermion, nine-orbital inequalities; supplementary row 52 equals the
printed three-fermion, ten-orbital row 82 restricted to an empty tenth orbital.
That identity is checked explicitly; it does not establish the historical
intent of the unprinted row.
