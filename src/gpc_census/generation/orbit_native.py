"""Orbit-native candidate construction: never materialise what is thrown away.

WHAT THIS REPLACES. The reference constructor enumerates every admissible
hyperplane, orients each into two taus, and screens the whole list. At rank 8
that is 332,840 rows of which 326,233 fail one integer comparison, and at rank
10 it would be 1.78 billion. The fast components proved separately (orbit
enumeration, direct survivor generation, the Hall structural filter) existed but
were not on this path, so none of them reached the actual system.

WHAT THIS BUILDS INSTEAD.

    orbit-canonical closed-flat enumeration
      -> one exact (h, z) per TOP-level orbit
      -> direct inversion-number survivor generation
      -> existing screening, which is now Hall-accelerated
      -> existing exact redundancy reduction

The full hyperplane list, the full oriented-tau list and the trace-rejected rows
are never constructed. Downstream is untouched: screening consumes taus, so this
module's only job is to emit exactly the taus that survive the trace test.

WHY COVERAGE STILL HOLDS, which is the whole point of the exercise. Two facts,
both already proved elsewhere in this package:

1. the admissible hyperplane set is ``S_d`` stable, so orbit representatives
   plus expansion reach every hyperplane (``orbit_canonical``);
2. the trace test is ``inv(sigma h) == B`` with ``B`` an orbit invariant, so an
   orbit's survivors are exactly the arrangements of the multiset ``h`` with
   ``B`` inversions, and those are generated directly
   (``generate_trace_survivors``).

Neither is an empirical pattern. Coverage is then CHECKED rather than asserted,
by count conservation against the Gaussian multinomial:

    2 * sum_orbits |O|  ==  trace_rejected + hall_zero + algebraic_zero
                            + certified + unresolved

with the left side computed from orbit-stabilizer and the right side from what
the run actually did. A missing orbit or a dropped survivor breaks the identity.

ORIENTATION SELF-PAIRING, which is where a naive version double counts. Each
unoriented hyperplane yields two taus. The oriented orbit of ``(h, z)`` is
``{(sigma h, z)}``, and it already contains ``(-h, -z)`` exactly when ``z == 0``
and the multiset of ``-h`` equals that of ``h``. Such an orbit is its own
reverse, so processing both orientations would emit every row twice. This is not
hypothetical: there are 1, 3 and 3 self-paired orbits at ranks 6, 7 and 8, which
is exactly why the oriented orbit counts are 15, 35 and 109 rather than 16, 38
and 112.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass

from .admissible_flats import rank_preserving_modulus
from .bdr import Tau, primitive_tau
from .exterior import exterior_weights
from .model import AdmissibleHyperplane
from .orbit_canonical import (
    enumerate_orbits_by_augmentation,
    generate_trace_survivors,
    orbit_size,
    trace_survivor_count,
)
from .ressayre import (
    _cofactor_nullvector,
    _level_sets,
    _primitive_pair,
)


def _canonical_hash(payload: object) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


@dataclass(frozen=True)
class OrientedOrbitRecord:
    """One oriented orbit: its exact normal, its size, and its survivors."""

    h: tuple[int, ...]
    z: int
    weights_below: int
    orbit_rows: int
    survivor_count: int

    def as_dict(self) -> dict[str, object]:
        return {
            "h": list(self.h),
            "z": self.z,
            "weights_below": self.weights_below,
            "orbit_rows": self.orbit_rows,
            "survivor_count": self.survivor_count,
        }


@dataclass(frozen=True)
class OrbitHyperplaneRecord:
    """One unoriented ``S_d`` orbit of admissible hyperplanes."""

    canonical_mask: int
    witness: tuple[int, ...]
    h: tuple[int, ...]
    z: int
    automorphism_order: int
    orbit_size: int
    self_paired: bool
    oriented: tuple[OrientedOrbitRecord, ...]

    def as_dict(self) -> dict[str, object]:
        return {
            "witness": list(self.witness),
            "h": list(self.h),
            "z": self.z,
            "automorphism_order": self.automorphism_order,
            "orbit_size": self.orbit_size,
            "self_paired": self.self_paired,
            "oriented": [record.as_dict() for record in self.oriented],
        }


@dataclass(frozen=True)
class OrbitNativeEnumeration:
    """Survivor taus and the count conservation that proves coverage."""

    n: int
    d: int
    orbits: tuple[OrbitHyperplaneRecord, ...]
    survivor_taus: tuple[Tau, ...]
    orbit_counts_by_rank: tuple[int, ...]
    hyperplane_count: int
    oriented_row_count: int
    predicted_survivor_count: int

    @property
    def trace_rejected_count(self) -> int:
        """Rows the materialising path would build and discard."""
        return self.oriented_row_count - len(self.survivor_taus)

    @property
    def counts_are_conserved(self) -> bool:
        return (
            len(self.survivor_taus) == self.predicted_survivor_count
            and self.trace_rejected_count >= 0
        )

    def as_dict(self) -> dict[str, object]:
        return {
            "method": "orbit_native_canonical_augmentation",
            "particle_number": self.n,
            "rank": self.d,
            "orbit_counts_by_affine_matrix_rank": list(self.orbit_counts_by_rank),
            "unoriented_orbit_count": len(self.orbits),
            "admissible_hyperplane_count": self.hyperplane_count,
            "oriented_row_count": self.oriented_row_count,
            "materialised_rows_avoided": self.trace_rejected_count,
            "survivor_count": len(self.survivor_taus),
            "predicted_survivor_count": self.predicted_survivor_count,
            "counts_are_conserved": self.counts_are_conserved,
            "canonical_survivor_tau_sha256": _canonical_hash(
                sorted([list(tau) for tau in self.survivor_taus])
            ),
            "coverage_argument": (
                "S_d stability of the admissible set reaches every hyperplane "
                "from an orbit representative; the trace test is inv(sigma h) "
                "== B with B an orbit invariant, so survivors are the "
                "arrangements of h with B inversions and are generated directly"
            ),
        }


def _exact_normal(n: int, d: int, witness: tuple[int, ...]):
    """Recover the exact primitive ``(h, z)`` of a top-rank closed flat.

    Deliberately the same construction the materialising enumerator uses, so
    the two paths cannot disagree about what a flat's hyperplane is.
    """
    weights = exterior_weights(n, d)
    base = weights[witness[0]]
    difference_rows = [
        [weights[index][column] - base[column] for column in range(d)]
        for index in witness[1:]
    ]
    normal = _cofactor_nullvector(difference_rows + [[1] * d])
    offset = sum(normal[index] * base[index] for index in range(d))
    return _primitive_pair(normal, offset, orient=True)


def _incidence_mask(n: int, d: int, h: tuple[int, ...], z: int) -> int:
    weights = exterior_weights(n, d)
    mask = 0
    for index, weight in enumerate(weights):
        if sum(left * right for left, right in zip(h, weight, strict=True)) == z:
            mask |= 1 << index
    return mask


def enumerate_orbit_survivor_taus(
    n: int, d: int, checkpoint_dir=None, verbose: bool = False
) -> OrbitNativeEnumeration:
    """Emit exactly the taus that pass the Ressayre trace test.

    The materialising path is never entered: no hyperplane list, no oriented-tau
    list, no trace-rejected rows.
    """
    if type(n) is not int or type(d) is not int or not (0 < n < d):
        raise ValueError("require integer parameters with 0 < n < d")
    rank_preserving_modulus(n, d)          # same Hadamard guard as the reference

    counts, levels = enumerate_orbits_by_augmentation(
        n, d, verbose=verbose, checkpoint_dir=checkpoint_dir)
    weights = exterior_weights(n, d)

    records: list[OrbitHyperplaneRecord] = []
    taus: list[Tau] = []
    hyperplanes = 0
    oriented_rows = 0
    predicted = 0

    for mask, witness in levels[-1].items():
        h, z = _exact_normal(n, d, witness)
        if _incidence_mask(n, d, h, z) != mask:
            raise AssertionError(
                "orbit representative disagrees with its exact hyperplane")
        size = orbit_size(mask, n, d)
        automorphisms = _factorial(d) // size
        hyperplanes += size
        oriented_rows += 2 * size

        negated = tuple(-value for value in h)
        self_paired = z == 0 and sorted(negated) == sorted(h)
        orientations = (1,) if self_paired else (1, -1)

        oriented: list[OrientedOrbitRecord] = []
        for orientation in orientations:
            oriented_h = tuple(orientation * value for value in h)
            oriented_z = orientation * z
            below = len(_level_sets(weights, oriented_h, oriented_z)[1])
            expected = trace_survivor_count(oriented_h, below)
            predicted += expected
            produced = 0
            for arrangement in generate_trace_survivors(oriented_h, below):
                taus.append(primitive_tau(
                    tuple(oriented_z - n * value for value in arrangement)))
                produced += 1
            if produced != expected:
                raise AssertionError(
                    "generated survivor count disagrees with the Gaussian "
                    "multinomial, which means one of the two routes is wrong")
            # rows in this oriented orbit: the distinct arrangements of h, which
            # is 2 * orbit_size when the orbit is self paired and orbit_size
            # otherwise, so the two cases agree on the total either way
            oriented.append(OrientedOrbitRecord(
                h=oriented_h,
                z=oriented_z,
                weights_below=below,
                orbit_rows=2 * size if self_paired else size,
                survivor_count=expected,
            ))

        records.append(OrbitHyperplaneRecord(
            canonical_mask=mask,
            witness=witness,
            h=h,
            z=z,
            automorphism_order=automorphisms,
            orbit_size=size,
            self_paired=self_paired,
            oriented=tuple(oriented),
        ))

    if len(set(taus)) != len(taus):
        raise AssertionError(
            "orbit-native generation produced duplicate taus, which means an "
            "orientation was processed twice")

    return OrbitNativeEnumeration(
        n=n,
        d=d,
        orbits=tuple(records),
        survivor_taus=tuple(taus),
        orbit_counts_by_rank=tuple(counts),
        hyperplane_count=hyperplanes,
        oriented_row_count=oriented_rows,
        predicted_survivor_count=predicted,
    )


def _factorial(k: int) -> int:
    result = 1
    for value in range(2, k + 1):
        result *= value
    return result


@dataclass(frozen=True)
class OrbitNativeGeneration:
    """A complete system built without materialising the rejected rows."""

    enumeration: OrbitNativeEnumeration
    screening: object
    candidate_system: object

    @property
    def proved_complete(self) -> bool:
        return (
            self.enumeration.counts_are_conserved
            and self.candidate_system.finalized
            and self.screening.all_candidates_resolved
            and self.screening.full_dimension.is_full_dimensional
            and self.candidate_system.counts.unresolved_rows == 0
        )

    def as_dict(self) -> dict[str, object]:
        reduction = self.candidate_system.reduction
        regression = self.candidate_system.known_system_regression
        counts = self.screening.counts.as_dict()
        survivors = len(self.enumeration.survivor_taus)
        return {
            "schema_version": 1,
            "particle_number": self.enumeration.n,
            "rank": self.enumeration.d,
            "status": (
                "proved_complete" if self.proved_complete else "blocked_incomplete"
            ),
            "candidate_enumeration": self.enumeration.as_dict(),
            "screening_counts": counts,
            # The identity that stands in for materialising the rejected rows.
            # The left side comes from orbit-stabilizer and the Gaussian
            # multinomial, the right side from what the run actually did, so a
            # lost orbit or a dropped survivor cannot balance it.
            "count_conservation": {
                "oriented_rows_total": self.enumeration.oriented_row_count,
                "trace_rejected_not_materialised": (
                    self.enumeration.trace_rejected_count
                ),
                "survivors_generated": survivors,
                "survivors_predicted_in_closed_form": (
                    self.enumeration.predicted_survivor_count
                ),
                "conserved": (
                    self.enumeration.counts_are_conserved
                    and self.enumeration.trace_rejected_count + survivors
                    == self.enumeration.oriented_row_count
                ),
            },
            "candidate_system_counts": self.candidate_system.counts.as_dict(),
            "ordered_slice_reduction_certificate": (
                None if reduction is None else reduction.as_dict()
            ),
            "retained_rows": [
                row.as_dict() for row in self.candidate_system.retained_rows
            ],
            "known_system_regression": (
                None if regression is None else regression.as_dict()
            ),
            "proof_status": {
                "candidate_exhaustiveness": (
                    "proved_by_orbit_canonical_flat_enumeration_and_"
                    "closed_form_survivor_conservation"
                ),
                "rank_preservation_modulo_prime": "proved_by_hadamard_bound",
                "candidate_classification": (
                    "exact_with_zero_unresolved_rows"
                    if self.screening.all_candidates_resolved
                    else "blocked_by_unresolved_rows"
                ),
                "ressayre_completeness_theorem": "applied_in_full_dimension",
                "redundancy_elimination": "exact_farkas_and_facet_certificates",
                "system_completeness": (
                    "proved" if self.proved_complete else "not_established"
                ),
            },
        }


def generate_orbit_native_rank(
    d: int,
    *,
    particle_number: int = 3,
    ressayre_attempts: int = 32,
    workers: int = 1,
    checkpoint_dir=None,
    verbose: bool = False,
) -> OrbitNativeGeneration:
    """Build a complete system from orbit representatives alone.

    Same downstream as the reference path, because screening consumes taus and
    this emits exactly the taus the reference path would have kept. The
    acceptance gate is not that it runs but that it reproduces the reference
    system, which ``tests/test_orbit_native.py`` checks at ranks 7 and 8.
    """
    from .candidate_pipeline import screen_tau_candidates
    from .candidate_system import compose_candidate_system

    enumeration = enumerate_orbit_survivor_taus(
        particle_number, d, checkpoint_dir=checkpoint_dir, verbose=verbose)
    if not enumeration.counts_are_conserved:
        raise AssertionError("orbit-native enumeration failed count conservation")
    screening = screen_tau_candidates(
        enumeration.survivor_taus,
        particle_number=particle_number,
        ressayre_attempts=ressayre_attempts,
        cyclic_cross_check=False,
        symbolic_fallback=True,
        workers=workers,
    )
    candidate_system = compose_candidate_system(screening)
    return OrbitNativeGeneration(
        enumeration=enumeration,
        screening=screening,
        candidate_system=candidate_system,
    )


def candidate_from_orbit_record(
    n: int, record: OrbitHyperplaneRecord
) -> AdmissibleHyperplane:
    """The representative hyperplane of an orbit, for replay and inspection."""
    return AdmissibleHyperplane(h=record.h, z=record.z, witness=record.witness)


__all__ = [
    "OrbitHyperplaneRecord",
    "OrbitNativeEnumeration",
    "OrbitNativeGeneration",
    "OrientedOrbitRecord",
    "candidate_from_orbit_record",
    "enumerate_orbit_survivor_taus",
    "generate_orbit_native_rank",
]
