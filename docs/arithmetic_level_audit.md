# Arithmetic-level audit of known GPC facet mechanisms

**Status:** post-hoc exploratory audit, not preregistered.

The audit combines the distinct mechanisms in the existing Levi-fusion pilot and held-out artifacts.

## Result

- Mechanisms tested: **83**
- Exact arithmetic progressions: **78**
- Arithmetic progressions with exactly one deleted point: **5**
- Mechanisms with more than one deleted point: **0**
- Cases where the progression step divides \(N\): **83/83**

Thus 78 of 83 mechanisms are exact arithmetic progressions, and the remaining five are one-point-deleted progressions.

## Elementary divisibility lemma

Suppose the distinct levels are \(a+gj\), with \(\gcd(a,g)=1\), and an \(N\)-particle zero-grade allocation exists. Then

\[
Na+gS=0
\]

for an integer \(S\). Therefore \(g\mid Na\), and primitiveness gives \(g\mid N\).

## Five exceptions

- (3,10): normalized levels `[0, 2, 3, 4, 5, 6, 7, 8]`, missing `[1]`.
- (4,10): normalized levels `[0, 1, 2, 3, 4, 5, 7]`, missing `[6]`.
- (4,9): normalized levels `[0, 1, 2, 3, 4, 5, 7]`, missing `[6]`.
- (5,10): normalized levels `[0, 1, 2, 3, 4, 5, 6, 7, 8, 10]`, missing `[9]`.
- (5,10): normalized levels `[0, 2, 3, 4, 5, 6, 7, 8, 9, 10]`, missing `[1]`.

## Constructor implication

A candidate mechanism can be parameterized by a divisor \(g\mid N\), an arithmetic interval (possibly with one hole), its multiplicity vector, and a zero-grade target. This replaces arbitrary integer normal multisets by a much smaller grammar, if the pattern survives an audit of all certified and rejected candidates.
