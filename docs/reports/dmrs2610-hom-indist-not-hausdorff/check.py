#!/usr/bin/env python3
"""Exact finite checks for the homomorphism distinguishing topology proposal.

Only Python's standard library is used. Finite checks do not prove the
unbounded sparse-moment lemma or the topological conclusion.
"""
from fractions import Fraction
from itertools import combinations


def graphs(n):
    edges = list(combinations(range(n), 2))
    for mask in range(1 << len(edges)):
        a = [[0] * n for _ in range(n)]
        for bit, (i, j) in enumerate(edges):
            if mask & (1 << bit):
                a[i][j] = a[j][i] = 1
        yield a


def cycle(n):
    assert n >= 3
    a = [[0] * n for _ in range(n)]
    for i in range(n):
        a[i][(i + 1) % n] = a[(i + 1) % n][i] = 1
    return a


def complete_bipartite(p, q, isolated=0):
    a = [[0] * (p + q + isolated) for _ in range(p + q + isolated)]
    for i in range(p):
        for j in range(p, p + q):
            a[i][j] = a[j][i] = 1
    return a


def disjoint_union(a, b):
    n, m = len(a), len(b)
    return [[a[i][j] if i < n and j < n else
             b[i - n][j - n] if i >= n and j >= n else 0
             for j in range(n + m)] for i in range(n + m)]


def connected(a):
    if not a:
        return False
    seen, todo = {0}, [0]
    while todo:
        i = todo.pop()
        for j in range(len(a)):
            if a[i][j] and j not in seen:
                seen.add(j)
                todo.append(j)
    return len(seen) == len(a)


def hom(guest, target):
    """Count all adjacency-preserving maps by exact backtracking, allowing repeats."""
    n, m = len(guest), len(target)
    earlier = [[j for j in range(i) if guest[i][j]] for i in range(n)]
    assigned = [0] * n

    def extend(i):
        if i == n:
            return 1
        total = 0
        for x in range(m):
            if all(target[x][assigned[j]] for j in earlier[i]):
                assigned[i] = x
                total += extend(i + 1)
        return total

    return extend(0)


def mul(a, b):
    n = len(a)
    return [[sum(a[i][k] * b[k][j] for k in range(n))
             for j in range(n)] for i in range(n)]


def traces(a, bound):
    n = len(a)
    power = [[int(i == j) for j in range(n)] for i in range(n)]
    ans = [n]
    for _ in range(bound):
        power = mul(power, a)
        ans.append(sum(power[i][i] for i in range(n)))
    return ans


def characteristic_polynomial(a):
    """Newton identities, returning coefficients in descending order."""
    t = traces(a, len(a))
    c = [1]
    for k in range(1, len(a) + 1):
        numerator = -sum(c[k - j] * t[j] for j in range(1, k + 1))
        assert numerator % k == 0
        c.append(numerator // k)
    return c


def determinant(a):
    m = [[Fraction(x) for x in row] for row in a]
    sign, result = 1, Fraction(1)
    for j in range(len(m)):
        pivot = next((i for i in range(j, len(m)) if m[i][j]), None)
        if pivot is None:
            return Fraction(0)
        if pivot != j:
            m[j], m[pivot] = m[pivot], m[j]
            sign *= -1
        z = m[j][j]
        result *= z
        for i in range(j + 1, len(m)):
            ratio = m[i][j] / z
            for k in range(j + 1, len(m)):
                m[i][k] -= ratio * m[j][k]
    return sign * result


def main():
    small = [a for n in range(5) for a in graphs(n)]
    guests = [a for a in small if connected(a)]
    targets = [a for a in small if len(a) <= 3]
    assert len(small) == 76 and len(guests) == 44 and len(targets) == 12
    assert hom([], []) == 1
    assert hom(cycle(4), []) == 0

    fidelity = 0
    # Source Section 3, equation labelled eq:connected (display (3.3)).
    for z in guests:
        for g in targets:
            for h in targets:
                assert hom(z, disjoint_union(g, h)) == hom(z, g) + hom(z, h)
                fidelity += 1
    print(f"FIDELITY: connected-guest disjoint-union identity: {fidelity} exact cases")
    print("DOMAIN: empty targets included; empty graph has one map to the empty graph")

    identities = 0
    for g in small:
        t = traces(g, 18)
        for m in range(3, 9):
            assert hom(cycle(m), g) == t[m]
            identities += 1
    print(f"CYCLE_TRACE: hom(C_m,G)=tr(A_G^m), 0<=|G|<=4, 3<=m<=8: {identities} cases")

    for p, q, r, s, iso, label in [
        (1, 4, 2, 2, 1, 'K1,4 versus C4 plus K1'),
        (1, 6, 2, 3, 2, 'K1,6 versus K2,3 plus two isolated vertices'),
    ]:
        g = complete_bipartite(p, q)
        h = complete_bipartite(r, s, iso)
        assert len(g) == len(h)
        assert sorted(map(sum, g)) != sorted(map(sum, h))
        expected_charpoly = [1, 0, -p * q] + [0] * (len(g) - 2)
        assert characteristic_polynomial(g) == characteristic_polynomial(h) == expected_charpoly
        for n in range(2, 6):
            expected = 2 * (p * q) ** n
            assert hom(cycle(2 * n), g) == hom(cycle(2 * n), h) == expected
        print(f"COSPECTRAL: {label}; charpoly={expected_charpoly}; even counts="
              f"{[2 * (p * q) ** n for n in range(2, 6)]}")

    # A different-order, noncospectral pair with the same positive squared spectrum.
    g, h = complete_bipartite(1, 4), complete_bipartite(2, 2)
    assert len(g) != len(h)
    assert characteristic_polynomial(g) != characteristic_polynomial(h)
    assert traces(g, 18)[1:] == traces(h, 18)[1:]
    print("ZERO_EIGENVALUES: K1,4 and C4 have different orders but equal positive-length cycle traces")

    # Finite evidence for the proposed sparse-moment implication; not an infinite-set test.
    signatures = {}
    for g in small:
        t = traces(g, 18)
        sparse = tuple(t[m] for m in [4, 10, 18])
        all_even = tuple(t[m] for m in range(4, 19, 2))
        if sparse in signatures:
            assert signatures[sparse] == all_even
        signatures[sparse] = all_even
    print(f"SPARSE_SMALL_CASES: 76 labelled targets; moments 4,10,18 determine moments 4,6,...,18; "
          f"{len(signatures)} signatures")

    vandermonde = 0
    for bases in [(1,), (1, 2), (1, 2, 5), (2, 3, 7), (1, 2, 4, 9)]:
        for step in [1, 2, 3]:
            exponents = [2 + step * i for i in range(len(bases))]
            assert determinant([[b ** n for b in bases] for n in exponents]) != 0
            vandermonde += 1
    assert determinant([[b ** n for b in (1, 2, 5)] for n in (2, 5, 9)]) != 0
    vandermonde += 1
    print(f"EXPONENTIAL_SUM: {vandermonde} exact positive-base nonsingular moment matrices")
    separators = 0
    for roots in [(), (2,), (3,), (2, 3)]:
        coefficients = [1]
        for n in roots:
            new = [0] * (len(coefficients) + 1)
            for j, c in enumerate(coefficients):
                new[j] -= (2 ** n) * c
                new[j + 1] += c
            coefficients = new
        positive, negative = [], []
        for j, c in enumerate(coefficients):
            for _ in range(abs(c)):
                star = complete_bipartite(1, 2 ** j)
                if c > 0:
                    positive = disjoint_union(positive, star)
                else:
                    negative = disjoint_union(negative, star)
        for n in range(2, 5):
            value = sum(c * (2 ** n) ** j for j, c in enumerate(coefficients))
            difference = hom(cycle(2 * n), positive) - hom(cycle(2 * n), negative)
            assert difference == 2 * value
            assert (difference == 0) == (n in roots)
            separators += 1
    print(f"FINITE_CYCLE_SEPARATORS: {separators} exact polynomial/star-union cases, including the empty set")
    print("LIMIT: these are finite arithmetic checks; the infinite-family theorem needs the proposed general proof")
    print("ALL_CHECKS_OK")


if __name__ == '__main__':
    main()
