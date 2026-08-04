# Stage 1 checkpoint: first-principles reproduction of the (3,7) system

## What is in this repository

The proof manifest `results/data/stage1_ressayre_3_7.json` and the independent
verifier `scripts/verify_stage1_ressayre_3_7.py`.

The generation package that produced the manifest is **not** in this
repository. The upstream write-up names
`src/gpc_census/generation/admissible_flats.py`,
`.../ressayre.py`, `.../reference_generator.py`,
`scripts/generate_fixed_n3_reference_rank.py` and
`scripts/verify_stage1_3_7.py`; none of those exist here, and the two modules
that were supplied import six siblings (`exterior`, `model`, `ressayre`,
`bdr`, `candidate_pipeline`, `candidate_system`) that were not. Until that
package lands, the manifest cannot be regenerated in this checkout, and the
regenerate-and-compare verifier described upstream cannot be written.

What can be done without it is to check the manifest against mathematics
rather than against its author, which is what the verifier here does.

## Result

The fixed-`N=3` rank-7 system reduces to four primitive inequalities:

\[
\begin{aligned}
-2\lambda_1+\lambda_2+\lambda_3+\lambda_4+\lambda_5
 -2\lambda_6-2\lambda_7 &\leq 0,\\
\lambda_1-2\lambda_2+\lambda_3+\lambda_4-2\lambda_5
 +\lambda_6-2\lambda_7 &\leq 0,\\
\lambda_1+\lambda_2-2\lambda_3-2\lambda_4+\lambda_5
 +\lambda_6-2\lambda_7 &\leq 0,\\
\lambda_1+\lambda_2-2\lambda_3+\lambda_4-2\lambda_5
 -2\lambda_6+\lambda_7 &\leq 0.
\end{aligned}
\]

These agree exactly, as an unordered set of oriented rows, with the
Altunbulak-Klyachko rank-7 table already stored in
`src/gpc_census/data/constraints.json`. That agreement is a regression check
on the manifest, not evidence of exhaustiveness.

## Exact affine-flat enumeration

The 35 weights of \(\bigwedge^3\mathbb C^7\), each augmented to
\((\omega,1)\), have rational rank 7, so an admissible affine hyperplane is a
closed flat of augmented rank 6. The generator builds closed flats one rank at
a time: if \(F\) is closed, two weights outside \(F\) give the same one-rank
extension exactly when their residuals modulo \(\operatorname{span}(F)\) are
projectively equal, and every rank-\(k+1\) closed flat \(G\) contains a
rank-\(k\) closed subflat \(F\) and some \(u\in G\setminus F\) with
\(G=\operatorname{cl}(F\cup\{u\})\).

| Augmented rank | Closed flats |
|---:|---:|
| 0 | 1 |
| 1 | 35 |
| 2 | 595 |
| 3 | 4,655 |
| 4 | 15,505 |
| 5 | 19,019 |
| 6 | 5,341 |

The verifier re-enumerates this lattice in exact rational arithmetic under
`--flats`, using no finite field at all, and reproduces every count.

## Why the finite-field closure is exact

Each augmented \(\bigwedge^n\) weight has \(n+1\) unit entries, so Hadamard
bounds every \(k\times k\) minor by \((n+1)^{k/2}\). The implementation
compares squares to avoid floating point: it requires \(p^2>(n+1)^d\). At
\(n=3\), \(d=7\) that is \(p^2>16384\), and \(257^2=66049\) satisfies it. A
rationally nonzero minor therefore has nonzero residue, so rank, span
membership, closure and projective residual equivalence over \(\mathbb F_p\)
are the characteristic-zero relations. This is a proof, not a probabilistic
check.

Note the doc-level shorthand that the prime must exceed \(2^7=128\) is the
same condition: \(\sqrt{(n+1)^d}=4^{3.5}=128\).

## Classification and reduction

Both orientations of all 5,341 hyperplanes give 10,682 candidate rows:

| Outcome | Rows |
|---|---:|
| Structural \(\lambda_7\geq0\) row | 1 |
| Exact trace-condition rejection | 10,133 |
| Exact zero tangent determinant | 499 |
| Nonstructural nonzero Ressayre determinant | 49 |
| Unresolved | 0 |

The 49 certified rows are reduced on the ordered trace-three Pauli slice:
45 are removed by exact Farkas certificates and 4 are retained with exact
rational facet witnesses.

## What the verifier checks, and what it cannot

Verified here, independently of the generator:

- the closed-flat lattice, re-enumerated exactly over \(\mathbb Q\);
- the Hadamard bound and the prime's sufficiency;
- the full-dimension certificate's seeded coefficients, regenerated from the
  recorded rule, and its skew and symmetric ranks against \(\binom d2\) and
  \(d(d+1)/2\);
- all 45 Farkas removals, replayed from their recorded multipliers, including
  that every multiplier is nonnegative and the residual slack is nonnegative;
- all 4 facet witnesses: each violates its own row by the recorded margin,
  satisfies every other retained row and the whole ordered slice, and has
  trace three;
- the relative interior point;
- the Ressayre on/below index splits, recomputed from \(h\) and \(z\) against
  the 35 weights, and \(\tau\) as the primitive form of \(z-N h\);
- the retained rows against the repository's stored table.

Not checkable without the generation package, and therefore not confirmed by
anything in this repository:

- the Ressayre tangent determinants recorded for each certified row;
- the 10,133 trace-condition rejections and the 499 symbolic zeros;
- the exhaustiveness of the candidate run, which is the load-bearing claim
  behind `"system_completeness": "proved"`;
- the cyclic-Schubert coefficient-one cross-checks.

The manifest's own `screening_counts` are worth reading carefully on the last
point: `cyclic_schubert_cross_checked_rows` is 0 and
`cyclic_schubert_deferred_rows` is 49, so the cross-check was deferred during
screening and appears only in the four retained rows afterwards.

## Reproduction

```sh
uv run python scripts/verify_stage1_ressayre_3_7.py           # certificate replay
uv run python scripts/verify_stage1_ressayre_3_7.py --flats   # plus exact lattice
```
