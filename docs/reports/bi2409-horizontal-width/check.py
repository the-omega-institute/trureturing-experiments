#!/usr/bin/env python3
"""Exact periodic-permutation checks for Bardakov--Iskra Question 6.2.

Standard-library integer arithmetic only. A common period N represents the
entire integer permutation by its values on 0,...,N-1; these generators fix
each N-block. Finite tests support, but do not replace, the general argument.
"""
from itertools import product
from math import lcm


def tau(m, r, s, x):
    assert 2 <= m and 0 <= r < s < m
    q, a = divmod(x, m)
    return q * m + (s if a == r else r if a == s else a)


def word_value(word, x):
    # Mathlib Perm multiplication: (f*g)(x) = f(g(x)).
    for m, r, s in reversed(word):
        x = tau(m, r, s, x)
    return x


def period(word):
    return lcm(*(m for m, _, _ in word)) if word else 1


def displacement(word):
    return {word_value(word, x) - x for x in range(period(word))}


def order(word):
    n = period(word)
    perm = [word_value(word, x) for x in range(n)]
    assert sorted(perm) == list(range(n))
    seen = set()
    answer = 1
    for start in range(n):
        if start in seen:
            continue
        x, length = start, 0
        while x not in seen:
            seen.add(x)
            length += 1
            x = perm[x]
        assert x == start
        answer = lcm(answer, length)
    return answer


def rho(n, x):
    assert n >= 2
    q, r = divmod(x, n)
    return q * n + n - 1 - r


def reversal_word(n):
    return tuple((n, r, n - 1 - r) for r in range(n // 2))


def ceil_log3(n):
    k, p = 0, 1
    while p < n:
        k, p = k + 1, p * 3
    return k


def main():
    # Fidelity: valid displayed orders in the proof of source Corollary 5.5,
    # plus the second example of Section 3 (orders are orientation invariant).
    examples = [
        (((2, 0, 1), (2, 0, 1)), 1),
        (((3, 0, 1), (3, 0, 2)), 3),
        (((2, 0, 1), (3, 0, 1)), 4),
        (((2, 0, 1), (3, 0, 2)), 6),
        (((3, 0, 1), (4, 0, 2)), 12),
        (((3, 0, 1), (4, 2, 3)), 6),
    ]
    for word, expected in examples:
        actual = order(word)
        assert actual == expected, (word, actual, expected)
        print(f"FIDELITY word={word} order={actual}")

    # The v1 displayed order-2 example has an unrelated typographical error.
    # Preserve its measured reading rather than silently treating it as valid.
    anomaly = ((2, 0, 1), (4, 0, 2))
    assert order(anomaly) == 4
    print("SOURCE_ANOMALY Corollary 5.5 displays order 2 for tau_0(2),1(2)*tau_0(4),2(4); exact order is 4")

    # Also reproduce the explicit cycle in the first Section-3 example.
    # Source acts on the right, so its written product reverses the
    # Mathlib-style composition order used elsewhere in this script.
    w = ((3, 0, 1), (2, 0, 1))
    assert [word_value(w, x) for x in range(6)] == [0, 1, 4, 2, 5, 3]
    for q in range(-5, 6):
        for r in range(6):
            assert word_value(w, 6*q+r) == 6*q + word_value(w, r)
    print("FIDELITY S6 cycle=(2,4,5,3); residues 0,1 fixed; negative blocks checked")

    generators = [(m, r, s) for m in range(2, 10)
                  for r in range(m) for s in range(r+1, m)]
    ds = {}
    for a in generators:
        m, r, s = a
        ds[a] = displacement((a,))
        assert ds[a] <= {0, s-r, r-s} and len(ds[a]) <= 3
        for x in range(-2*m, 2*m+1):
            assert tau(*a, tau(*a, x)) == x
            assert tau(*a, x+m) == tau(*a, x)+m
    pair_orders = set()
    for a, b in product(generators, repeat=2):
        d = displacement((a, b))
        sums = {u+v for u in ds[a] for v in ds[b]}
        assert d <= sums and len(d) <= len(ds[a])*len(ds[b])
        pair_orders.add(order((a, b)))
    assert pair_orders == {1, 2, 3, 4, 6, 12}
    print(f"PAIRS generators={len(generators)} moduli=2..9 ordered_pairs={len(generators)**2}")
    print(f"PAIRS orders={sorted(pair_orders)}; displacement inclusion and cardinal bound hold")

    small = [a for a in generators if a[0] <= 4]
    count = 0
    for k in range(5):
        for word in product(small, repeat=k):
            assert len(displacement(word)) <= 3**k
            count += 1
    assert displacement(()) == {0}
    print(f"WORDS moduli=2..4 lengths=0..4 count={count}; |D|<=3^length including empty word")

    for n in range(2, 129):
        word = reversal_word(n)
        expected = {n-1-2*r for r in range(n)}
        assert len(expected) == n
        assert {rho(n, r)-r for r in range(n)} == expected
        for q in range(-3, 4):
            for r in range(n):
                x = q*n+r
                y = rho(n, x)
                assert rho(n, y) == x
                assert word_value(word, x) == y
                assert (y == x) == (n % 2 == 1 and r == n//2)
        assert rho(n, 0) != 0
    print("FAMILY n=2..128 q=-3..3 all residues: factorization, involution, exact D, fixed points checked")

    for b in range(8):
        n = 3**b+1
        word = reversal_word(n)
        d = {rho(n, r)-r for r in range(n)}
        assert n % 2 == 0 and len(d) == n and 0 not in d
        for r in range(n):
            assert word_value(word, r) == rho(n, r)
        assert len(d) > 3**b and ceil_log3(n) > b
        assert all(len(d) > 3**k for k in range(b+1))
        print(f"WITNESS B={b} n={n} displacements={len(d)} explicit_factors={len(word)} lower_bound={ceil_log3(n)} fixed_point_free=True")
    assert {rho(4, r)-r for r in range(4)} == {-3, -1, 1, 3}
    print("ALL_CHECKS_OK")


if __name__ == '__main__':
    main()
