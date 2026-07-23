"""Murnaghan-Nakayama character tables for S_n (exact integers). Verified by
dimension sum-of-squares = n! and row orthogonality."""
import itertools
from functools import lru_cache
from math import prod, factorial


def partitions(n, mx=None):
    if mx is None:
        mx = n
    if n == 0:
        yield ()
        return
    for k in range(min(n, mx), 0, -1):
        for rest in partitions(n - k, k):
            yield (k,) + rest


def cycle_types(n):
    return list(partitions(n))


def class_size(mu, n):
    from collections import Counter
    c = Counter(mu)
    z = prod(factorial(c[k]) * (k ** c[k]) for k in c)
    return factorial(n) // z


@lru_cache(maxsize=None)
def mn_char(lam, mu):
    """chi^lambda(mu) via Murnaghan-Nakayama recursion. lam,mu sorted-desc tuples."""
    if sum(lam) == 0:
        return 1
    k = mu[0]
    rest = mu[1:]
    total = 0
    # remove all border strips of size k from lam
    for strip in border_strips(lam, k):
        newlam, height = strip
        total += (-1) ** height * mn_char(newlam, rest)
    return total


def border_strips(lam, k):
    """yield (lambda_minus_strip, height) for every size-k border strip of lam."""
    lam = list(lam)
    r = len(lam)
    # represent by rows; a border strip removal: pick a connected skew rim of size k.
    # Use the standard beta-number / rim-hook method.
    # beta_i = lam_i + (r-1-i)
    beta = [lam[i] + (r - 1 - i) for i in range(r)]
    bset = set(beta)
    for i in range(r):
        b = beta[i]
        if b - k >= 0 and (b - k) not in bset:
            # remove hook: beta[i] -> beta[i]-k
            newbeta = beta[:]; newbeta[i] = b - k
            # height = number of beta_j strictly between b-k and b
            height = sum(1 for j in range(r) if b - k < beta[j] < b)
            nb = sorted(newbeta, reverse=True)
            newlam = tuple(nb[j] - (r - 1 - j) for j in range(r))
            newlam = tuple(x for x in newlam if x > 0)
            yield (newlam, height)


def char_table(n):
    lams = cycle_types(n)   # irreps <-> partitions
    mus = cycle_types(n)    # classes  <-> partitions
    tab = {}
    for lam in lams:
        for mu in mus:
            tab[(lam, mu)] = mn_char(lam, mu)
    return lams, mus, tab


def verify(n):
    lams, mus, tab = char_table(n)
    dims = [tab[(lam, (1,) * n)] for lam in lams]
    ok_dim = sum(d * d for d in dims) == factorial(n)
    # row orthogonality: sum_mu |class| chi_a(mu) chi_b(mu) = n! delta_ab
    ok_orth = True
    for a in lams:
        for b in lams:
            s = sum(class_size(mu, n) * tab[(a, mu)] * tab[(b, mu)] for mu in mus)
            expect = factorial(n) if a == b else 0
            if s != expect:
                ok_orth = False
    return dims, ok_dim, ok_orth


if __name__ == "__main__":
    for n in (5, 6):
        dims, ok_dim, ok_orth = verify(n)
        print(f"S{n}: irrep dims {sorted(dims)}  sum sq = {sum(d*d for d in dims)} (= {factorial(n)}? {ok_dim})  "
              f"orthogonality {ok_orth}")
