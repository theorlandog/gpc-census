"""First-principles generalized Pauli constraint generation."""

from .admissible_flats import (
    FlatHyperplaneEnumeration,
    enumerate_admissible_hyperplanes_by_flats,
    rank_preserving_modulus,
)

from .candidate_pipeline import (
    CandidateScreeningCounts,
    CandidateScreeningResult,
    ScreenedCandidateRow,
    StructuralPauliRow,
    SymbolicFallbackOutcome,
    UnresolvedCandidatesError,
    screen_bdr_export,
    screen_tau_candidates,
)
from .candidate_system import (
    BLOCKED_KNOWN_MISMATCH_STATUS,
    BLOCKED_UNRESOLVED_STATUS,
    FINALIZED_KNOWN_STATUS,
    FINALIZED_RELATIVE_STATUS,
    CandidateSystemCounts,
    CandidateSystemInvariantError,
    CandidateSystemResult,
    FinalCandidateRow,
    KnownSystemRegression,
    compose_candidate_system,
)
from .full_dimension import (
    FullDimensionCertificate,
    full_dimension_certificate,
    require_full_dimension,
    verify_full_dimension_certificate,
)
from .model import GeneratedSystem, IntegralConstraint, RessayreCertificate
from .redundancy import (
    CertificateSearchError,
    FacetWitnessCertificate,
    FarkasRemovalCertificate,
    OrderedSliceReductionCertificate,
    certify_ordered_slice_redundancy,
    ordered_pauli_inequalities,
    verify_ordered_slice_reduction_certificate,
)
from .reference_generator import (
    FixedN3ReferenceGeneration,
    generate_fixed_n3_reference_rank,
    verify_reference_manifest,
)
from .system import generate_3_6

__all__ = [
    "BLOCKED_KNOWN_MISMATCH_STATUS",
    "BLOCKED_UNRESOLVED_STATUS",
    "FINALIZED_KNOWN_STATUS",
    "FINALIZED_RELATIVE_STATUS",
    "CandidateScreeningCounts",
    "CandidateScreeningResult",
    "CandidateSystemCounts",
    "CandidateSystemInvariantError",
    "CandidateSystemResult",
    "CertificateSearchError",
    "FacetWitnessCertificate",
    "FlatHyperplaneEnumeration",
    "FarkasRemovalCertificate",
    "FinalCandidateRow",
    "FixedN3ReferenceGeneration",
    "FullDimensionCertificate",
    "GeneratedSystem",
    "IntegralConstraint",
    "KnownSystemRegression",
    "OrderedSliceReductionCertificate",
    "RessayreCertificate",
    "ScreenedCandidateRow",
    "StructuralPauliRow",
    "SymbolicFallbackOutcome",
    "UnresolvedCandidatesError",
    "certify_ordered_slice_redundancy",
    "compose_candidate_system",
    "enumerate_admissible_hyperplanes_by_flats",
    "full_dimension_certificate",
    "generate_3_6",
    "generate_fixed_n3_reference_rank",
    "ordered_pauli_inequalities",
    "rank_preserving_modulus",
    "require_full_dimension",
    "screen_bdr_export",
    "screen_tau_candidates",
    "verify_full_dimension_certificate",
    "verify_ordered_slice_reduction_certificate",
    "verify_reference_manifest",
]
