#!/usr/bin/env python3
"""Exact finite checks for the repeated-singleton MAO inequality.

Only integer counts and fractions.Fraction carry mathematical evidence.
The unrestricted family argument is in the accompanying preregistration;
the finite checks here do not certify an infinite theorem or a Lean proof.
"""

from collections import Counter
from fractions import Fraction
from itertools import combinations, product
from math import comb, prod


def falling(x, r):
    """The source product, including its zero factor when r > x >= 0."""
    assert x >= 0 and r >= 0
    return prod(x - j for j in range(r))


def patterns(T, t):
    return tuple(sum(1 << i for i in S) for S in combinations(range(T), t))


def weight(n, sizes, slots):
    r = len(slots)
    result = 1
    for i, m in enumerate(sizes):
        k = sum((S >> i) & 1 for S in slots)
        result *= falling(m, k) * falling(n - m, r - k)
    return result


def transversal(n, sizes, levels):
    T = len(sizes)
    choices = [tuple(S for t in B for S in patterns(T, t)) for B in levels]
    return sum(weight(n, sizes, slots) for slots in product(*choices))


def moment(n, sizes, levels):
    r = len(levels)
    assert 1 <= r <= n and len(sizes) >= 2
    assert all(0 <= m <= n for m in sizes)
    denominator = falling(n, r) ** (len(sizes) - 1)
    assert denominator > 0
    return Fraction(transversal(n, sizes, levels), denominator)


def exact_moment(n, sizes, r, t):
    return moment(n, sizes, ((t,),) * r)


def enumerate_counts(n, sizes):
    """Independent original sample space: actual prescribed-size subsets."""
    T = len(sizes)
    samples = [patterns(n, m) for m in sizes]
    exact = [Counter() for _ in range(T + 1)]
    threshold = [Counter() for _ in range(T + 1)]
    configurations = 0
    for chosen_sets in product(*samples):
        occupancy = [sum((S >> a) & 1 for S in chosen_sets) for a in range(n)]
        histogram = Counter(occupancy)
        cumulative = 0
        for t in range(T, -1, -1):
            count = histogram[t]
            cumulative += count
            exact[t][count] += 1
            threshold[t][cumulative] += 1
        configurations += 1
    assert configurations == prod(comb(n, m) for m in sizes)
    return configurations, exact, threshold


def expectation(hist, denominator, f):
    return Fraction(sum(count * f(x) for x, count in hist.items()), denominator)


def raw_from_factorial(mu):
    return (mu[0], mu[0] + mu[1], mu[0] + 3 * mu[1] + mu[2],
            mu[0] + 7 * mu[1] + 6 * mu[2] + mu[3])


def source_fidelity():
    # arXiv:2608.23998v1 section 8 explicitly states this exhaustive example.
    n, sizes = 9, (3, 2, 2, 1)
    total, exact, threshold = enumerate_counts(n, sizes)
    assert total == 979776
    print(f"FIDELITY source section 8: n={n}, sizes={sizes}, configurations={total}")
    checked = 0
    for t in range(1, 5):
        for label, hist, B in (("exact", exact[t], (t,)),
                               ("threshold", threshold[t], tuple(range(t, 5)))):
            factorial = tuple(moment(n, sizes, (B,) * r) for r in range(1, 5))
            for r, value in enumerate(factorial, 1):
                assert value == expectation(hist, total, lambda x: falling(x, r))
                checked += 1
            raw = raw_from_factorial(factorial)
            for r, value in enumerate(raw, 1):
                assert value == expectation(hist, total, lambda x: x ** r)
                checked += 1
            variance = raw[1] - raw[0] ** 2
            assert variance == expectation(hist, total, lambda x: x ** 2) - raw[0] ** 2
            checked += 1
            print(f"  t={t} {label}: mean={raw[0]}, variance={variance}, raw1_to_4={raw}")
    for r in range(1, n + 1):
        # Theorem 6 / equation (92), the actually established fully occupied case.
        assert exact_moment(n, sizes, r, 4) <= exact_moment(n, sizes, 1, 4) ** r
    print(f"FIDELITY_OK: {checked} exact identities; fully occupied inequality r=1,...,{n}")


def smallest_witness():
    n, sizes, t = 2, (1, 1, 1, 1), 2
    S = patterns(4, t)
    assert len(S) == 6
    surviving = []
    for A, B in product(S, repeat=2):
        w = weight(n, sizes, (A, B))
        assert w == int(B == (A ^ 15))
        if w:
            surviving.append((A, B))
    assert len(surviving) == 6
    a = exact_moment(n, sizes, 1, t)
    b = exact_moment(n, sizes, 2, t)
    total, hist, _ = enumerate_counts(n, sizes)
    assert total == 16 and hist[t] == Counter({0: 10, 2: 6})
    assert a == b == Fraction(3, 4)
    assert a * a == Fraction(9, 16) < b
    variance = expectation(hist[t], total, lambda x: x ** 2) - a * a
    assert variance == Fraction(15, 16) > a
    print("WITNESS n=r=2, T=4, all m_i=1, t=2")
    print(f"  ordered_pattern_pairs=36, surviving_pairs={len(surviving)}, denominators=8,8")
    print(f"  E[X]={a}, E[(X)_2]={b}, E[X]^2={a*a}, excess={b-a*a}, variance={variance}")


def family_checks():
    # Two independent finite counts: all central pattern pairs, then all samples.
    for k in range(1, 6):
        T, sizes = 2 * k, (1,) * (2 * k)
        u = Fraction(comb(T, k), 2 ** (T - 1))
        assert exact_moment(2, sizes, 1, k) == u
        assert exact_moment(2, sizes, 2, k) == u
        all_bits = (1 << T) - 1
        central = patterns(T, k)
        for A, B in product(central, repeat=2):
            assert weight(2, sizes, (A, B)) == int(B == (A ^ all_bits))
    for k in range(1, 8):
        T = 2 * k
        total, exact, _ = enumerate_counts(2, (1,) * T)
        assert exact[k] == Counter({2: comb(T, k), 0: 2 ** T - comb(T, k)})
        u = Fraction(comb(T, k), 2 ** (T - 1))
        assert expectation(exact[k], total, lambda x: x) == u
        assert expectation(exact[k], total, lambda x: falling(x, 2)) == u
        variance = expectation(exact[k], total, lambda x: x ** 2) - u ** 2
        assert variance == 2 * u - u ** 2
        if k >= 2:
            assert 0 < u < 1 and u ** 2 < u and variance > u
        print(f"  FAMILY k={k}: u_k={u}, u_k-u_k^2={u-u*u}, variance={variance}")
    assert Fraction(comb(2, 1), 2) == 1
    for k in range(1, 101):
        u = Fraction(comb(2 * k, k), 2 ** (2 * k - 1))
        v = Fraction(comb(2 * k + 2, k + 1), 2 ** (2 * k + 1))
        assert v * (2 * k + 2) == u * (2 * k + 1)
        if k >= 2:
            assert 0 < u < 1 and u ** 2 < u
    print("FAMILY_OK: pattern pairs k=1,...,5; all samples k=1,...,7; recurrence k=1,...,100")
    print("INFINITE_SCOPE: finite checks support the outlined bijection/induction; no infinite proof claimed")


def main():
    source_fidelity()
    smallest_witness()
    family_checks()
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
