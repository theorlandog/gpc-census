"""Proof-carrying fixed-``N=3`` reference generation by affine flats."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass

from .admissible_flats import (
    FlatHyperplaneEnumeration,
    enumerate_admissible_hyperplanes_by_flats,
)
from .bdr import Tau, primitive_tau
from .candidate_pipeline import CandidateScreeningResult, screen_tau_candidates
from .candidate_system import CandidateSystemResult, compose_candidate_system
from .model import AdmissibleHyperplane


def _canonical_hash(payload: object) -> str:
    encoded = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def tau_from_ressayre_candidate(
    candidate: AdmissibleHyperplane,
    particle_number: int,
    orientation: int,
) -> Tau:
    """Convert one oriented ``H.lambda >= z`` pair to ``tau.lambda <= 0``."""
    if orientation not in (-1, 1):
        raise ValueError("orientation must be -1 or 1")
    h = tuple(orientation * value for value in candidate.h)
    z = orientation * candidate.z
    return primitive_tau(tuple(z - particle_number * value for value in h))


def oriented_taus(
    enumeration: FlatHyperplaneEnumeration,
) -> tuple[Tau, ...]:
    """Return both primitive orientations of every enumerated hyperplane."""
    rows = tuple(
        tau_from_ressayre_candidate(candidate, enumeration.n, orientation)
        for candidate in enumeration.hyperplanes
        for orientation in (1, -1)
    )
    if len(set(rows)) != len(rows):
        raise AssertionError("oriented hyperplane conversion produced duplicate taus")
    return rows


def _status_partition(
    screening: CandidateScreeningResult,
) -> dict[str, dict[str, object]]:
    grouped: dict[str, list[list[int]]] = {}
    for row in screening.candidates:
        grouped.setdefault(row.status, []).append(list(row.tau))
    grouped["structural_pauli_lower_bound"] = [
        list(row.tau) for row in screening.structural_rows
    ]
    return {
        status: {
            "count": len(rows),
            "canonical_tau_sha256": _canonical_hash(sorted(rows)),
        }
        for status, rows in sorted(grouped.items())
    }


@dataclass(frozen=True)
class FixedN3ReferenceGeneration:
    """A complete exhaustive run and its exact downstream certificates."""

    enumeration: FlatHyperplaneEnumeration
    taus: tuple[Tau, ...]
    screening: CandidateScreeningResult
    candidate_system: CandidateSystemResult

    @property
    def proved_complete(self) -> bool:
        return (
            self.candidate_system.finalized
            and self.screening.all_candidates_resolved
            and self.screening.full_dimension.is_full_dimensional
            and self.candidate_system.counts.unresolved_rows == 0
        )

    def as_dict(self) -> dict[str, object]:
        """Return a compact manifest whose verifier regenerates the run."""
        reduction = self.candidate_system.reduction
        regression = self.candidate_system.known_system_regression
        return {
            "schema_version": 1,
            "particle_number": 3,
            "rank": self.enumeration.d,
            "status": (
                "proved_complete" if self.proved_complete else "blocked_incomplete"
            ),
            "candidate_enumeration": {
                **self.enumeration.as_dict(),
                "method": "exact_closed_affine_flat_lattice",
                "orientation_count": len(self.taus),
                "canonical_oriented_tau_sha256": _canonical_hash(
                    sorted([list(tau) for tau in self.taus])
                ),
                "coverage_argument": (
                    "every rank-k closed flat is generated from a rank-(k-1) "
                    "closed subflat and one projective residual class"
                ),
            },
            "screening_counts": self.screening.counts.as_dict(),
            "candidate_status_partition": _status_partition(self.screening),
            "ambient_full_dimension_certificate": (
                self.screening.full_dimension.as_dict()
            ),
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
                "candidate_exhaustiveness": "proved_by_exact_flat_enumeration",
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


def generate_fixed_n3_reference_rank(
    d: int,
    *,
    ressayre_attempts: int = 32,
    workers: int = 1,
    checkpoint_dir: str | None = None,
) -> FixedN3ReferenceGeneration:
    """Generate a complete fixed-``N=3`` system from affine flats.

    This is a coverage-first reference path. It is practical at ranks 7 and 8
    and is intended as the exact oracle against which faster pruning backends
    are proved. Resource exhaustion at a larger rank has no mathematical
    meaning.

    ``checkpoint_dir`` makes the flat enumeration resumable across runs, which
    is what a rank whose cost is hours rather than gigabytes needs.
    """
    if type(d) is not int or not 7 <= d <= 13:
        raise ValueError(
            "the reference composition currently supports ranks 7 through 13"
        )
    enumeration = enumerate_admissible_hyperplanes_by_flats(3, d, checkpoint_dir)
    taus = oriented_taus(enumeration)
    screening = screen_tau_candidates(
        taus,
        ressayre_attempts=ressayre_attempts,
        cyclic_cross_check=False,
        symbolic_fallback=True,
        workers=workers,
    )
    candidate_system = compose_candidate_system(screening)
    result = FixedN3ReferenceGeneration(
        enumeration=enumeration,
        taus=taus,
        screening=screening,
        candidate_system=candidate_system,
    )
    if not result.proved_complete:
        raise RuntimeError(
            "reference generation did not close every exact completeness gate"
        )
    return result


def verify_reference_manifest(
    payload: dict[str, object],
    *,
    workers: int = 1,
) -> bool:
    """Regenerate and compare one compact reference manifest exactly."""
    try:
        if payload.get("schema_version") != 1:
            return False
        rank = payload.get("rank")
        if type(rank) is not int:
            return False
        regenerated = generate_fixed_n3_reference_rank(rank, workers=workers)
        return regenerated.as_dict() == payload
    except (AssertionError, RuntimeError, TypeError, ValueError):
        return False


__all__ = [
    "FixedN3ReferenceGeneration",
    "generate_fixed_n3_reference_rank",
    "oriented_taus",
    "tau_from_ressayre_candidate",
    "verify_reference_manifest",
]
