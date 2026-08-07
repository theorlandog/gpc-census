"""Guards for the parabolic reformulation and its refuted speedup.

The write-up rests on a NEGATIVE result, so the tests have to keep it negative:
a later change that made the incremental traversal look like a win would convert
a recorded refutation into a confirmation. The positive half, the parabolic
identification itself, is pinned separately because it is what survives.
"""
import collections
import itertools

import pytest

from gpc_census.generation import orbit_canonical as oc


def _inversions(sequence):
    return sum(1 for i in range(len(sequence))
               for j in range(i + 1, len(sequence))
               if sequence[j] < sequence[i])


# --------------------------------------------------------------------------
# the identification, which is correct
# --------------------------------------------------------------------------

@pytest.mark.parametrize("h", [(0, 1, 1, 1, 2, 2), (0, 0, 1, 1, 2), (0, 1, 2, 3)])
def test_arrangements_are_the_parabolic_quotient(h):
    d = len(h)
    young = [g for g in itertools.permutations(range(d))
             if all(h[g[i]] == h[i] for i in range(d))]
    order = 1
    for multiplicity in collections.Counter(h).values():
        for value in range(2, multiplicity + 1):
            order *= value
    assert len(young) == order
    arrangements = set(itertools.permutations(h))
    factorial = 1
    for value in range(2, d + 1):
        factorial *= value
    assert len(arrangements) == factorial // order


@pytest.mark.parametrize("h", [(0, 1, 1, 1, 2, 2), (0, 0, 1, 1, 2)])
def test_minimal_coset_length_is_the_inversion_number_only_when_h_increases(h):
    """The convention is load bearing and the other direction is FALSE.

    Sorted decreasingly, the base arrangement is the LONGEST representative, not
    the shortest, and neither the minimum nor the maximum coset length matches.
    """
    def cosets(base):
        grouped = collections.defaultdict(list)
        for g in itertools.permutations(range(len(base))):
            grouped[tuple(base[g[i]] for i in range(len(base)))].append(g)
        return grouped

    increasing = tuple(sorted(h))
    grouped = cosets(increasing)
    assert all(min(_inversions(g) for g in members) == _inversions(arrangement)
               for arrangement, members in grouped.items())

    decreasing = tuple(sorted(h, reverse=True))
    grouped = cosets(decreasing)
    assert not all(min(_inversions(g) for g in members) == _inversions(arrangement)
                   for arrangement, members in grouped.items())


def test_the_gaussian_multinomial_is_the_poincare_polynomial():
    """Convention free, and already what trace_survivor_count reads."""
    for h in [(0, 1, 1, 1, 2, 2), (0, 0, 1, 1, 2), (3, 1, 2)]:
        distribution = collections.Counter(
            _inversions(a) for a in set(itertools.permutations(h)))
        polynomial = oc.inversion_generating_function(h)
        assert list(polynomial) == [distribution.get(i, 0)
                                    for i in range(len(polynomial))]


# --------------------------------------------------------------------------
# the obstruction, which must stay an obstruction
# --------------------------------------------------------------------------

def test_no_two_equal_length_arrangements_are_adjacent():
    """A fixed-length Gray code cannot exist.

    An adjacent transposition changes Coxeter length by exactly one, so two
    distinct length-B representatives are never one swap apart. If this ever
    found a pair, the impossibility argument in the write-up would be wrong.
    """
    examined = preserving = 0
    for h in [(0, 1, 1, 2), (0, 0, 1, 1, 2), (0, 1, 1, 1, 2, 2), (0, 1, 2, 3, 3)]:
        arrangements = set(itertools.permutations(h))
        for arrangement in arrangements:
            for position in range(len(h) - 1):
                swapped = list(arrangement)
                swapped[position], swapped[position + 1] = (
                    swapped[position + 1], swapped[position])
                swapped = tuple(swapped)
                if swapped in arrangements and swapped != arrangement:
                    examined += 1
                    if _inversions(swapped) == _inversions(arrangement):
                        preserving += 1
    assert examined > 500, examined
    assert preserving == 0
