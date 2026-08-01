---
title: |
  Supplementary Material for “Symbolic pure-state preimages for tabulated vertex spectra of fermionic one-body marginals through ten orbitals”
---

# Scope of the supplementary files

The accompanying machine-readable archive contains the finite objects used in the manuscript's four main computational claims:

1. symbolic one-body reduced density matrix identities for 799 pure states;
2. finite no-root-distinct certificates for 12 spectra;
3. rational reconstruction of the vertices of nine supplied half-space descriptions; and
4. an ordered concordance between the rank-9 and rank-10 coefficient tables and Appendix A of Altunbulak's dissertation.

The archive does not reproduce the exploratory CP-SAT and CBC searches behind the remaining 144 solver-negative classification labels. Those labels are metadata in the full research repository and enter no theorem, proposition, or corollary in the manuscript. Fiber geometry, Hamiltonian dynamics, and rank-11 exploratory calculations are also outside the archive.

# Data objects

## Symbolic state records

The file `data/states.jsonl` has one JSON object for each row of the nine vertex tables. Its top-level fields are:

| field | meaning |
|:---|:---|
| `system` | the pair $(N,d)$ |
| `index` | zero-based row number in the corresponding vertex table |
| `classified` | the construction-search class retained as dataset metadata |
| `integer_form` | the canonical integer numerator of the target spectrum |
| `denominator` | its common denominator |
| `closed_form` | determinant support and symbolic amplitudes |

Orbital indices in the JSON arrays start at zero, whereas formulas in the manuscript start at one. The `closed_form.support_dets` array lists the occupied orbitals of each determinant. The corresponding `pretty` expressions are parsed as exact SymPy expressions; floating discovery amplitudes, timings, solver statuses, and internal transport paths are not included in this Paper 1 projection.

The state verifier applies $a_i^\dagger a_j$ directly to the determinant arrays, including fermionic insertion and removal signs. It proves normalization and compares the characteristic polynomial of the reconstructed 1-RDM with
$$
\prod_{k=1}^{d}\left(x-\frac{n_k}{D}\right)
$$
coefficient by coefficient. The bundle checker independently matches every `(system, index, integer_form, denominator)` key to one row of `data/vertices/`.

## No-root-distinct certificates

The file `data/interference_certificates.json` contains 12 finite proofs: the 11 negative cases at $(3,8)$ and the spectrum $v_B$ at $(4,9)$. Each record stores the target, integral Farkas separators, the equal-eigenvalue symmetry blocks, and a branch directed acyclic graph.

For a separator $y$, the checker recomputes the positive clause
$$
P_y=\{T:y\mathbin{\cdot}(\mathbf 1_T,1)>0\}
$$
and verifies $y\mathbin{\cdot}(v,1)>0$. At each branch node it checks one-hop independence of the selected determinants, confirms that the cited clause is still unhit, expands orbit representatives under the recorded symmetry group, and proves complete coverage of every compatible choice. Dead leaves have no compatible clause member. The checker also verifies acyclicity and reachability from the root.

## Half-space and vertex tables

The file `data/constraints.json` contains the nine normalized rational H-descriptions used as input. `data/vertices/` contains exactly nine corresponding rational V-tables. The H-to-V verifier begins with the ordered Pauli simplex, imposes the tabulated half-spaces in their stored order, and performs every intersection, rank, feasibility, and equality test with `fractions.Fraction`. It then compares the reconstructed vertices with the supplied tables in both directions.

This calculation proves vertex exhaustion for the supplied H-descriptions. It does not by itself prove representation-theoretic validity of the rank-10 inequalities, which remains the explicit hypothesis in the manuscript's conditional ten-orbital corollary.

## Source concordance

The file `data/constraint_source_map.json` identifies the dissertation PDF by SHA-256 and maps every printed coefficient row to a PDF page, thesis page, and row position. The checker compares the extraction with both `data/source_rows.json` and `data/constraints.json`.

All 490 printed rank-9 and rank-10 rows match: 51 at $(3,9)$, 60 at $(4,9)$, and 379 across the three rank-10 systems. Appendix A states that $(3,9)$ has 52 inequalities but prints 51. The final supplied $(3,9)$ row equals printed $(3,10)$ row 82 after setting $\lambda_{10}=0$. The archive verifies this algebraic identity but does not claim that it establishes the historical intent of the omitted row.

# Reproduction

The archive records its source commit and tool versions in `SOURCE.json`. `SHA256SUMS` covers every other file. From the extracted archive, create the pinned environment outside the archive directory so that the checksum namespace remains unchanged:

```sh
export UV_PROJECT_ENVIRONMENT=../gpc-paper1-venv
uv sync --locked --no-install-project
uv run --no-sync python scripts/verify_paper1_bundle_standalone.py .
uv run --no-sync python scripts/verify_states_standalone.py data/states.jsonl
uv run --no-sync python scripts/verify_interference_certificates_standalone.py \
  data/interference_certificates.json
uv run --no-sync python scripts/verify_vertex_exhaustion_standalone.py \
  --constraints data/constraints.json --vertices data/vertices
```

The dissertation PDF is not redistributed. After obtaining the bitstream identified in the source map, reproduce and compare the entire page-level concordance with:

```sh
uv run --no-sync python scripts/verify_constraint_source_standalone.py \
  dissertation.pdf \
  --source-rows data/source_rows.json \
  --normalized data/constraints.json \
  --expected-map data/constraint_source_map.json
```

The verification programs reject unknown files in the archive, missing or incomplete checksums, relabeled state or certificate targets, altered source rows, and disagreement between the recorded page map and a fresh extraction.

# Archive manifest

The Paper 1 archive consists of the following groups:

- `README.md`, `SOURCE.json`, `LICENSE`, `.python-version`, `pyproject.toml`, and `uv.lock`;
- `SHA256SUMS`;
- the five top-level files under `data/` described above and the nine explicit files under `data/vertices/`; and
- five standalone programs under `scripts/` for state reconstruction, no-design proof checking, H-to-V reconstruction, source concordance, and cross-link verification.

No network access is needed for the symbolic state, certificate, H-to-V, or bundle checks after the environment has been created. The source-concordance check additionally requires the separately obtained dissertation PDF and a local `pdftotext` executable.
