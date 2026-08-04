"""Exact cyclic-Schubert coefficients for exterior-power inequalities.

This module implements the coefficient test needed by Klyachko's
one-particle marginal inequalities for ``Exterior^n(C^d)``.  For an affine
candidate

``sum_i coefficients[i] * lambda[i] <= rhs``

let ``a`` be the coefficients in weakly decreasing order, and let ``v`` be
the minimal permutation satisfying

``sum_i a[i] * lambda[v[i]] = sum_i coefficients[i] * lambda[i]``.

The essential cyclic target has length

``m = #{T in choose([d], n) : sum(a[i] for i in T) > rhs}``.

When ``length(v) == m`` and ``rhs`` is an exterior-weight value, its pulled
back Schubert polynomial is

``P(x) = product_T sum(i in T, x[i])``

over the same upper subsets.  The desired Schubert coefficient is the
degree-zero integer ``partial_v P``.

The exact evaluator uses no symbolic polynomial expansion.  If ``F_k`` is
the result after the first ``k`` divided differences in the fixed reduced
word, then at a permuted integral evaluation point ``b[p]`` it uses

``F_k(b[p]) = (F_(k-1)(b[p]) - F_(k-1)(s_i b[p]))
               / (b[p[i]] - b[p[i + 1]])``.

Every ``F_k`` lies in ``Z[x]`` because type-A divided differences preserve
integral polynomials.  Each division is therefore exact.  Once the degree
has fallen to zero, the final evaluation is the exact integer coefficient,
not a probabilistic test or a residue reconstruction.

The modular evaluator is a faster independent certificate.  A nonzero
residue modulo its recorded prime proves characteristic-zero nonvanishing.
It is not used to decide whether the integer coefficient equals one.

For ``S`` memoized evaluation states, ``m`` upper subsets, and particle
number ``n``, evaluation takes ``O(S + Lmn)`` arithmetic operations, where
``L`` is the number of leaf states, and ``O(S)`` memory.  In the worst case
``S`` is exponential in the reduced-word length.  Memoization identifies
all repeated permutation states and is substantially smaller on the known
fixed-``n`` systems.
"""

from __future__ import annotations

import itertools
from collections import defaultdict
from dataclasses import dataclass
from functools import lru_cache
from typing import Sequence


DEFAULT_MODULUS = 2_147_483_647


Permutation = tuple[int, ...]
Subset = tuple[int, ...]


@dataclass(frozen=True)
class CyclicSchubertProblem:
    """Deterministic combinatorial data for one cyclic-Schubert test.

    Permutations, divided-difference indices, and exterior subsets are all
    zero-based.  ``source_permutation`` is the minimal representative through
    equal-coefficient blocks.  This makes the problem data independent of
    incidental sorting choices.
    """

    particle_number: int
    coefficients: tuple[int, ...]
    rhs: int
    dominant_coefficients: tuple[int, ...]
    source_permutation: Permutation
    divided_difference_word: tuple[int, ...]
    upper_weight_subsets: tuple[Subset, ...]
    threshold_is_weight: bool

    @property
    def degree_matches(self) -> bool:
        """Whether the source and cyclic target have the same length."""
        return len(self.divided_difference_word) == len(self.upper_weight_subsets)

    @property
    def is_well_formed(self) -> bool:
        """Whether the cyclic coefficient is a degree-zero intersection number."""
        return self.degree_matches and self.threshold_is_weight

    def as_dict(self) -> dict[str, object]:
        """Return a deterministic, JSON-serializable problem payload."""
        return {
            "particle_number": self.particle_number,
            "coefficients": list(self.coefficients),
            "rhs": self.rhs,
            "dominant_coefficients": list(self.dominant_coefficients),
            "source_permutation_zero_based": list(self.source_permutation),
            "divided_difference_word_zero_based": list(self.divided_difference_word),
            "upper_weight_subsets_zero_based": [
                list(subset) for subset in self.upper_weight_subsets
            ],
            "threshold_is_weight": self.threshold_is_weight,
            "degree_matches": self.degree_matches,
        }


@dataclass(frozen=True)
class ExactCyclicSchubertCertificate:
    """Replayable exact-integer coefficient certificate."""

    problem: CyclicSchubertProblem
    evaluation_base: tuple[int, ...]
    coefficient: int
    evaluation_states: int

    @property
    def coefficient_one(self) -> bool:
        """Return the exact coefficient-one acceptance decision."""
        return self.coefficient == 1

    def as_dict(self) -> dict[str, object]:
        """Return a deterministic, JSON-serializable certificate payload."""
        return {
            "kind": "exact_integer_cyclic_schubert",
            "problem": self.problem.as_dict(),
            "evaluation_base": list(self.evaluation_base),
            "coefficient": self.coefficient,
            "coefficient_one": self.coefficient_one,
            "evaluation_states": self.evaluation_states,
        }


@dataclass(frozen=True)
class ModularCyclicSchubertCertificate:
    """Replayable finite-field certificate for the same coefficient."""

    problem: CyclicSchubertProblem
    modulus: int
    evaluation_base: tuple[int, ...]
    residue: int
    evaluation_states: int

    @property
    def proves_nonzero(self) -> bool:
        """Whether this residue proves characteristic-zero nonvanishing."""
        return self.residue != 0

    def as_dict(self) -> dict[str, object]:
        """Return a deterministic, JSON-serializable certificate payload."""
        return {
            "kind": "modular_cyclic_schubert",
            "problem": self.problem.as_dict(),
            "modulus": self.modulus,
            "evaluation_base_modulus": list(self.evaluation_base),
            "residue": self.residue,
            "proves_nonzero": self.proves_nonzero,
            "evaluation_states": self.evaluation_states,
        }


@dataclass(frozen=True)
class CyclicSchubertVerdict:
    """Exact acceptance verdict, with an optional modular cross-check."""

    problem: CyclicSchubertProblem
    status: str
    exact_certificate: ExactCyclicSchubertCertificate | None
    modular_certificate: ModularCyclicSchubertCertificate | None

    @property
    def coefficient(self) -> int | None:
        """Return the exact integer coefficient when the problem is well formed."""
        if self.exact_certificate is None:
            return None
        return self.exact_certificate.coefficient

    @property
    def coefficient_one(self) -> bool:
        """Return true exactly when the characteristic-zero coefficient is one."""
        return self.coefficient == 1

    def as_dict(self) -> dict[str, object]:
        """Return a deterministic, JSON-serializable verdict payload."""
        return {
            "status": self.status,
            "coefficient": self.coefficient,
            "coefficient_one": self.coefficient_one,
            "problem": self.problem.as_dict(),
            "exact_certificate": (
                None if self.exact_certificate is None else self.exact_certificate.as_dict()
            ),
            "modular_certificate": (
                None if self.modular_certificate is None else self.modular_certificate.as_dict()
            ),
        }


def _validate_candidate(
    coefficients: Sequence[int], rhs: int, particle_number: int
) -> tuple[int, ...]:
    if not isinstance(particle_number, int) or isinstance(particle_number, bool):
        raise TypeError("particle_number must be an integer")
    if particle_number <= 0:
        raise ValueError("particle_number must be positive")
    if not isinstance(rhs, int) or isinstance(rhs, bool):
        raise TypeError("rhs must be an integer")
    values = tuple(coefficients)
    if len(values) < particle_number:
        raise ValueError("the coefficient count must be at least particle_number")
    if any(not isinstance(value, int) or isinstance(value, bool) for value in values):
        raise TypeError("all coefficients must be integers")
    return values


def _minimal_assignment_permutation(
    dominant: tuple[int, ...], coefficients: tuple[int, ...]
) -> Permutation:
    positions: dict[int, list[int]] = defaultdict(list)
    for index, value in enumerate(coefficients):
        positions[value].append(index)
    used: dict[int, int] = defaultdict(int)
    permutation: list[int] = []
    for value in dominant:
        occurrence = used[value]
        permutation.append(positions[value][occurrence])
        used[value] = occurrence + 1
    return tuple(permutation)


def _inversion_count(permutation: Permutation) -> int:
    return sum(
        permutation[left] > permutation[right]
        for left in range(len(permutation))
        for right in range(left + 1, len(permutation))
    )


def _right_reduction_word(permutation: Permutation) -> tuple[int, ...]:
    """Return ``i`` with ``v s_i[0] ... s_i[-1] = identity``.

    The first adjacent inversion is removed at every step.  This fixes one
    deterministic reduced word and makes the serialized certificate stable.
    """
    current = list(permutation)
    word: list[int] = []
    while True:
        for index in range(len(current) - 1):
            if current[index] > current[index + 1]:
                current[index], current[index + 1] = current[index + 1], current[index]
                word.append(index)
                break
        else:
            break
    if current != list(range(len(permutation))):
        raise AssertionError("source data is not a permutation")
    if len(word) != _inversion_count(permutation):
        raise AssertionError("the adjacent-swap word is not reduced")
    return tuple(word)


def prepare_cyclic_schubert_problem(
    coefficients: Sequence[int], rhs: int, *, particle_number: int = 3
) -> CyclicSchubertProblem:
    """Construct the exact combinatorial problem for one candidate row.

    The construction is invariant under positive integral rescaling and under
    adding a multiple of the fixed-trace equation to the affine row.
    """
    values = _validate_candidate(coefficients, rhs, particle_number)
    dominant = tuple(sorted(values, reverse=True))
    source = _minimal_assignment_permutation(dominant, values)
    word = _right_reduction_word(source)
    exterior_subsets = tuple(
        itertools.combinations(range(len(values)), particle_number)
    )
    subset_values = tuple(
        sum(dominant[index] for index in subset) for subset in exterior_subsets
    )
    upper = tuple(
        subset
        for subset, value in zip(exterior_subsets, subset_values, strict=True)
        if value > rhs
    )
    return CyclicSchubertProblem(
        particle_number=particle_number,
        coefficients=values,
        rhs=rhs,
        dominant_coefficients=dominant,
        source_permutation=source,
        divided_difference_word=word,
        upper_weight_subsets=upper,
        threshold_is_weight=rhs in subset_values,
    )


def _validate_problem(problem: CyclicSchubertProblem) -> None:
    reconstructed = prepare_cyclic_schubert_problem(
        problem.coefficients,
        problem.rhs,
        particle_number=problem.particle_number,
    )
    if problem != reconstructed:
        raise ValueError("problem data is inconsistent with its candidate row")
    if not problem.degree_matches:
        raise ValueError("source length does not match the cyclic target length")
    if not problem.threshold_is_weight:
        raise ValueError("rhs is not an exterior-weight value")


def _evaluation_base(size: int, supplied: Sequence[int] | None) -> tuple[int, ...]:
    base = tuple(range(size)) if supplied is None else tuple(supplied)
    if len(base) != size:
        raise ValueError("evaluation base has the wrong length")
    if any(not isinstance(value, int) or isinstance(value, bool) for value in base):
        raise TypeError("evaluation base values must be integers")
    if len(set(base)) != size:
        raise ValueError("evaluation base values must be pairwise distinct")
    return base


def exact_cyclic_schubert_certificate(
    problem: CyclicSchubertProblem,
    *,
    evaluation_base: Sequence[int] | None = None,
) -> ExactCyclicSchubertCertificate:
    """Compute and certify the exact integer ``partial_v P``.

    Correctness follows by induction on the divided-difference word.  The
    memoized value at ``(k, p)`` is the integral polynomial after ``k``
    divided differences evaluated at the point whose coordinate ``i`` is
    ``base[p[i]]``.  The recurrence is precisely the definition of the
    type-A divided difference.  Exact divisibility is checked at every state.
    Degree matching makes the final polynomial constant, so its value is the
    desired integer coefficient at every pairwise-distinct base point.
    """
    _validate_problem(problem)
    size = len(problem.coefficients)
    base = _evaluation_base(size, evaluation_base)
    word = problem.divided_difference_word
    upper = problem.upper_weight_subsets

    @lru_cache(maxsize=None)
    def evaluate(depth: int, permutation: Permutation) -> int:
        if depth == 0:
            value = 1
            for subset in upper:
                value *= sum(base[permutation[index]] for index in subset)
            return value
        index = word[depth - 1]
        swapped = list(permutation)
        swapped[index], swapped[index + 1] = swapped[index + 1], swapped[index]
        numerator = evaluate(depth - 1, permutation) - evaluate(depth - 1, tuple(swapped))
        denominator = base[permutation[index]] - base[permutation[index + 1]]
        quotient, remainder = divmod(numerator, denominator)
        if remainder != 0:
            raise ArithmeticError("integral divided difference did not divide exactly")
        return quotient

    identity = tuple(range(size))
    coefficient = evaluate(len(word), identity)
    return ExactCyclicSchubertCertificate(
        problem=problem,
        evaluation_base=base,
        coefficient=coefficient,
        evaluation_states=evaluate.cache_info().currsize,
    )


def _is_prime_below_2_64(value: int) -> bool:
    """Deterministically test primality for an unsigned 64-bit integer."""
    if value < 2:
        return False
    small_primes = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    for prime in small_primes:
        if value % prime == 0:
            return value == prime
    odd_part = value - 1
    power_of_two = 0
    while odd_part % 2 == 0:
        power_of_two += 1
        odd_part //= 2
    for base in (2, 325, 9_375, 28_178, 450_775, 9_780_504, 1_795_265_022):
        if base % value == 0:
            continue
        witness = pow(base, odd_part, value)
        if witness in (1, value - 1):
            continue
        for _ in range(power_of_two - 1):
            witness = witness * witness % value
            if witness == value - 1:
                break
        else:
            return False
    return True


def modular_cyclic_schubert_certificate(
    problem: CyclicSchubertProblem,
    *,
    modulus: int = DEFAULT_MODULUS,
    evaluation_base: Sequence[int] | None = None,
) -> ModularCyclicSchubertCertificate:
    """Compute a deterministic finite-field nonvanishing certificate.

    Reduction modulo a prime commutes with every integral divided difference.
    Pairwise-distinct base residues make all displayed denominators invertible.
    The returned residue is therefore the exact integer coefficient modulo the
    recorded prime.  A nonzero residue rigorously proves that coefficient is
    nonzero over the integers, but residue one alone does not prove that the
    integer coefficient equals one.
    """
    _validate_problem(problem)
    if not isinstance(modulus, int) or isinstance(modulus, bool):
        raise TypeError("modulus must be an integer")
    if modulus >= 2**64 or not _is_prime_below_2_64(modulus):
        raise ValueError("modulus must be a prime smaller than 2^64")
    size = len(problem.coefficients)
    integral_base = _evaluation_base(size, evaluation_base)
    base = tuple(value % modulus for value in integral_base)
    if len(set(base)) != size:
        raise ValueError("evaluation base residues must be pairwise distinct")
    word = problem.divided_difference_word
    upper = problem.upper_weight_subsets
    inverse_difference = {
        (left, right): pow((base[left] - base[right]) % modulus, modulus - 2, modulus)
        for left in range(size)
        for right in range(size)
        if left != right
    }

    @lru_cache(maxsize=None)
    def evaluate(depth: int, permutation: Permutation) -> int:
        if depth == 0:
            value = 1
            for subset in upper:
                linear_value = sum(base[permutation[index]] for index in subset)
                value = value * linear_value % modulus
            return value
        index = word[depth - 1]
        swapped = list(permutation)
        swapped[index], swapped[index + 1] = swapped[index + 1], swapped[index]
        difference = evaluate(depth - 1, permutation) - evaluate(depth - 1, tuple(swapped))
        return (
            difference * inverse_difference[permutation[index], permutation[index + 1]]
        ) % modulus

    identity = tuple(range(size))
    residue = evaluate(len(word), identity)
    return ModularCyclicSchubertCertificate(
        problem=problem,
        modulus=modulus,
        evaluation_base=base,
        residue=residue,
        evaluation_states=evaluate.cache_info().currsize,
    )


def certify_cyclic_schubert(
    coefficients: Sequence[int],
    rhs: int,
    *,
    particle_number: int = 3,
    include_modular: bool = True,
    modulus: int = DEFAULT_MODULUS,
) -> CyclicSchubertVerdict:
    """Return an exact coefficient-one verdict for one candidate inequality.

    Degree or threshold failures are ordinary negative verdicts.  A well-formed
    problem is accepted only when its exact integer coefficient equals one.
    The modular calculation, enabled by default, is an independent replayable
    cross-check and never substitutes for exact acceptance.
    """
    problem = prepare_cyclic_schubert_problem(
        coefficients,
        rhs,
        particle_number=particle_number,
    )
    if not problem.degree_matches:
        return CyclicSchubertVerdict(problem, "degree_mismatch", None, None)
    if not problem.threshold_is_weight:
        return CyclicSchubertVerdict(problem, "rhs_not_exterior_weight", None, None)
    exact = exact_cyclic_schubert_certificate(problem)
    modular = (
        modular_cyclic_schubert_certificate(problem, modulus=modulus)
        if include_modular
        else None
    )
    if modular is not None and exact.coefficient % modulus != modular.residue:
        raise AssertionError("exact and modular cyclic-Schubert evaluations disagree")
    status = "coefficient_one" if exact.coefficient_one else "coefficient_nonunit"
    return CyclicSchubertVerdict(problem, status, exact, modular)


def verify_exact_cyclic_schubert_certificate(
    certificate: ExactCyclicSchubertCertificate,
) -> bool:
    """Replay an exact certificate from its serialized mathematical data."""
    try:
        return certificate == exact_cyclic_schubert_certificate(
            certificate.problem,
            evaluation_base=certificate.evaluation_base,
        )
    except (ArithmeticError, TypeError, ValueError):
        return False


def verify_modular_cyclic_schubert_certificate(
    certificate: ModularCyclicSchubertCertificate,
) -> bool:
    """Replay a modular certificate from its modulus, base, and problem data."""
    try:
        return certificate == modular_cyclic_schubert_certificate(
            certificate.problem,
            modulus=certificate.modulus,
            evaluation_base=certificate.evaluation_base,
        )
    except (TypeError, ValueError):
        return False


__all__ = [
    "DEFAULT_MODULUS",
    "CyclicSchubertProblem",
    "CyclicSchubertVerdict",
    "ExactCyclicSchubertCertificate",
    "ModularCyclicSchubertCertificate",
    "certify_cyclic_schubert",
    "exact_cyclic_schubert_certificate",
    "modular_cyclic_schubert_certificate",
    "prepare_cyclic_schubert_problem",
    "verify_exact_cyclic_schubert_certificate",
    "verify_modular_cyclic_schubert_certificate",
]
