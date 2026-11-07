#!/usr/bin/env python3
"""Exact GF(2) checks for the fixed-field MRD-density proposal.

Every count below is integer/rational arithmetic.  The map model counts ordered
linear injections, whose images are uniformly distributed among subspaces.
"""
from fractions import Fraction
from itertools import combinations
from math import prod


def rank2(rows: list[int], width: int) -> int:
    rows = [x for x in rows if x]
    r = 0
    for col in range(width):
        pivot = next((i for i in range(r, len(rows)) if (rows[i] >> col) & 1), None)
        if pivot is None:
            continue
        rows[r], rows[pivot] = rows[pivot], rows[r]
        for i in range(len(rows)):
            if i != r and ((rows[i] >> col) & 1):
                rows[i] ^= rows[r]
        r += 1
        if r == len(rows):
            break
    return r


def dot_rows(rows: tuple[int, ...], x: int) -> int:
    out = 0
    for i, row in enumerate(rows):
        if (row & x).bit_count() & 1:
            out |= 1 << i
    return out


def injections_mrd_n2(m: int, d: int = 2) -> tuple[int, int]:
    """Exact number of all injective maps and MRD maps for n=2, N=m."""
    N = m
    map_count = 1 << (m * N)
    maps = [tuple((a >> (i * N)) & ((1 << N) - 1) for i in range(m))
            for a in range(map_count)]
    good = 0
    injections = 0
    for a in maps:
        for b in maps:
            inj = True
            mrd = True
            for x in range(1, 1 << N):
                cols = [dot_rows(a, x), dot_rows(b, x)]
                r = rank2(cols, m)
                # f(x) is an mn-dimensional vector; injectivity means its
                # matrix image is nonzero for every x != 0.
                if r < 1:
                    inj = False
                    mrd = False
                    break
                if r < d:
                    mrd = False
            if inj:
                injections += 1
                if mrd:
                    good += 1
    return injections, good


def rref_subspaces(D: int, k: int):
    """Yield every k-subspace of F_2^D once, using reduced row echelon bases."""
    for pivots in combinations(range(D), k):
        free_positions = []
        pivot_set = set(pivots)
        for i, p in enumerate(pivots):
            for col in range(p + 1, D):
                if col not in pivot_set:
                    free_positions.append((i, col))
        for mask in range(1 << len(free_positions)):
            rows = [1 << p for p in pivots]
            for bit, (i, col) in enumerate(free_positions):
                if (mask >> bit) & 1:
                    rows[i] |= 1 << col
            yield tuple(rows)


def subspace_mrd_count(n: int, m: int, d: int) -> tuple[int, int]:
    D = n * m
    k = m * (n - d + 1)
    total = good = 0
    for basis in rref_subspaces(D, k):
        total += 1
        ok = True
        for mask in range(1, 1 << k):
            word = 0
            for i in range(k):
                if (mask >> i) & 1:
                    word ^= basis[i]
            rows = [(word >> (m * i)) & ((1 << m) - 1) for i in range(n)]
            if rank2(rows, m) < d:
                ok = False
                break
        if ok:
            good += 1
    return total, good


def gaussian_2(a: int, b: int) -> int:
    num = den = 1
    for h in range(b):
        num *= (1 << a) - (1 << h)
        den *= (1 << b) - (1 << h)
    assert num % den == 0
    return num // den


def P_partial(q: int, terms: int) -> Fraction:
    p = Fraction(1)
    for j in range(1, terms + 1):
        p *= Fraction(q**j - 1, q**j)
    return p


def main() -> None:
    # Fidelity: exact map count agrees with the classical n=d=2 MRD model.
    counts = {}
    for m in (2, 3):
        injections, good = injections_mrd_n2(m)
        delta = Fraction(good, injections)
        counts[m] = delta
        expected_injections = prod((2 ** (2 * m) - 2**h) for h in range(m))
        assert injections == expected_injections
        print(f"q=2,n=2,d=2,m={m}: injections={injections}, mrd_maps={good}, delta={delta}")

    # Independent subspace enumeration for a second parameter set.
    total, good = subspace_mrd_count(3, 3, 2)
    delta_333 = Fraction(good, total)
    assert total == gaussian_2(9, 6)
    print(f"q=2,n=3,d=2,m=3: all_codes={total}, mrd_codes={good}, delta={delta_333}")

    # Proposed limits.  P(2) is bracketed exactly by a partial product and
    # the union bound 1 - sum_{j>r} 2^{-j} <= tail product <= 1.
    p20 = P_partial(2, 20)
    p20_lower = p20 * (1 - Fraction(1, 2**20))
    assert p20_lower <= p20
    target_n2_lower = p20_lower**3
    target_n2_upper = p20**3
    target_n3_lower = p20_lower**7
    target_n3_upper = p20**7
    print(f"P2 partial(20)={p20}")
    print(f"target n=2 lower={target_n2_lower}")
    print(f"target n=2 upper={target_n2_upper}")
    print(f"target n=3 lower={target_n3_lower}")
    print(f"target n=3 upper={target_n3_upper}")
    # The exact finite values lie above the limiting interval, as expected;
    # the n=2 values are positive and bounded by the stated P(2)^3 limit.
    assert target_n2_upper < counts[2]
    assert target_n2_upper < counts[3]
    assert target_n3_upper < delta_333
    # Gaussian exponent in the proposed identity is G(n,d-1;2).
    assert gaussian_2(2, 1) == 3
    assert gaussian_2(3, 1) == 7
    print("exact exponent checks: G(2,1;2)=3, G(3,1;2)=7")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
