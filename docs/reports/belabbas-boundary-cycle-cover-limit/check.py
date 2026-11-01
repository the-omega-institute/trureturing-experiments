#!/usr/bin/env python3
"""Exact finite fidelity and boundary-family checks; no asymptotic proof is claimed."""
from fractions import Fraction as Q
from functools import lru_cache
from itertools import combinations, permutations, product
from math import factorial, prod


def compositions(n, q):
    if q == 1:
        yield (n,)
    else:
        for first in range(n + 1):
            for rest in compositions(n - first, q - 1):
                yield (first,) + rest


def strongly_connected(q, edges):
    for start in range(q):
        seen = {start}
        while True:
            larger = seen | {j for i, j in edges if i in seen}
            if larger == seen:
                break
            seen = larger
        if len(seen) != q:
            return False
    return True


def neighbors(q, edges, subset):
    return {j for i, j in edges if i in subset}


def hall(q, edges, m):
    for mask in range(1, (1 << q) - 1):
        subset = {i for i in range(q) if mask >> i & 1}
        if sum(m[i] for i in subset) > sum(m[j] for j in neighbors(q, edges, subset)):
            return False
    return True


def perfect_matching(q, edges, m):
    """Independent bitmask matching search in the complete directed blow-up's lift."""
    types = tuple(i for i in range(q) for _ in range(m[i]))
    k = len(types)
    rows = tuple(sum(1 << b for b in range(k) if (types[a], types[b]) in edges)
                 for a in range(k))

    @lru_cache(None)
    def visit(used):
        a = used.bit_count()
        if a == k:
            return True
        choices = rows[a] & ~used
        while choices:
            bit = choices & -choices
            if visit(used | bit):
                return True
            choices -= bit
        return False
    return visit(0)


def source_fidelity():
    """Lemma 7 and Proposition 1: Hall feasibility iff a complete blow-up has a cover."""
    skeletons = cases = 0
    for q in (2, 3):
        arcs = tuple(product(range(q), repeat=2))
        for bits in range(1 << len(arcs)):
            edges = frozenset(e for j, e in enumerate(arcs) if bits >> j & 1)
            if not strongly_connected(q, edges):
                continue
            skeletons += 1
            for n in range(1, 7):
                for m in compositions(n, q):
                    assert hall(q, edges, m) == perfect_matching(q, edges, m), (q, edges, m)
                    cases += 1
    print(f'SOURCE_FIDELITY: {skeletons} strongly connected skeletons; {cases} count vectors; n=1..6')


def family_edges(q):
    return frozenset({(0, 0)} | {(i, i + 1) for i in range(q - 1)} |
                     {(i, 0) for i in range(1, q)})


def simple_cycles(q, edges):
    """Canonicalize by starting at the smallest vertex; orientation is preserved."""
    result = set()
    for length in range(1, q + 1):
        for support in combinations(range(q), length):
            for tail in permutations(support[1:]):
                cycle = (support[0],) + tail
                if all((cycle[j], cycle[(j + 1) % length]) in edges for j in range(length)):
                    result.add(cycle)
    return result


def rank_rational(rows):
    a = [[Q(x) for x in row] for row in rows]
    rank = 0
    for col in range(len(a[0])):
        pivot = next((r for r in range(rank, len(a)) if a[r][col]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        pivot_value = a[rank][col]
        a[rank] = [x / pivot_value for x in a[rank]]
        for r in range(len(a)):
            if r != rank:
                c = a[r][col]
                a[r] = [x - c * y for x, y in zip(a[r], a[rank])]
        rank += 1
    return rank


def coefficients(q, edges, mask):
    subset = {i for i in range(q) if mask >> i & 1}
    out = neighbors(q, edges, subset)
    return tuple(int(i in out) - int(i in subset) for i in range(q))


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def family_geometry():
    for q in range(2, 7):
        edges = family_edges(q)
        assert strongly_connected(q, edges)
        cycles = simple_cycles(q, edges)
        assert cycles == {tuple(range(j)) for j in range(1, q + 1)}
        vertices = [tuple(Q(int(i < j), j) for i in range(q)) for j in range(1, q + 1)]
        assert rank_rational([[v[i] - vertices[0][i] for i in range(q)]
                              for v in vertices[1:]]) == q - 1
        # Exact inverse barycentric formula: lambda_j = j*(x_j-x_{j+1}), x_{q+1}=0.
        for n in range(1, 9):
            for m in compositions(n, q):
                x = tuple(Q(a, n) for a in m)
                lambdas = tuple((j + 1) * (x[j] - (x[j + 1] if j + 1 < q else 0))
                                for j in range(q))
                assert sum(lambdas) == 1
                assert all(sum(lambdas[j] * vertices[j][i] for j in range(q)) == x[i]
                           for i in range(q))
                ordered = all(m[i] >= m[i + 1] for i in range(q - 1))
                assert (min(lambdas) >= 0) == ordered == hall(q, edges, m)
        xstar = (Q(1, q),) * q
        covariance = tuple(tuple(Q(int(i == j), q) - Q(1, q * q)
                                 for j in range(q)) for i in range(q))
        assert all(sum(row) == 0 for row in covariance)
        assert all(covariance[i][j] == covariance[(i + 1) % q][(j + 1) % q]
                   for i in range(q) for j in range(q))
        active = [coefficients(q, edges, mask) for mask in range(1, (1 << q) - 1)
                  if dot(coefficients(q, edges, mask), xstar) == 0]
        assert active
        for alpha in active:
            variance = sum(alpha[i] * covariance[i][j] * alpha[j]
                           for i in range(q) for j in range(q))
            assert variance == sum(xstar[i] * alpha[i] ** 2 for i in range(q)) > 0
        for i in range(q):
            for j in range(i):
                assert covariance[i][i] + covariance[j][j] - 2 * covariance[i][j] == Q(2, q)
        accepted = total = 0
        for order in permutations(range(q)):
            z = tuple(Q(2 * a - (q - 1), 2) for a in order)
            cone = all(dot(alpha, z) >= 0 for alpha in active)
            assert cone == all(z[i] >= z[i + 1] for i in range(q - 1))
            accepted += cone
            total += 1
        assert accepted == 1 and total == factorial(q)
        print(f'FAMILY_GEOMETRY: q={q}; cycles={q}; dimension={q-1}; active={len(active)}; chamber_fraction={Q(accepted,total)}')


def attainable_sizes(q, m):
    """Pack prefix cycles independently of the Hall inequalities."""
    @lru_cache(None)
    def visit(remaining):
        result = {0}
        for j in range(1, q + 1):
            if all(remaining[i] > 0 for i in range(j)):
                smaller = tuple(a - int(i < j) for i, a in enumerate(remaining))
                result.update(j + k for k in visit(smaller))
        return frozenset(result)
    return visit(tuple(m))


def all_sizes_and_probability():
    cases = 0
    for q in range(2, 6):
        edges = family_edges(q)
        for n in range(1, 11):
            total_weight = ordered_weight = distinct_weight = ordered_distinct_weight = 0
            for m in compositions(n, q):
                denominator = prod(factorial(a) for a in m)
                assert factorial(n) % denominator == 0
                weight = factorial(n) // denominator
                sizes = attainable_sizes(q, m)
                all_sizes = all(k in sizes for k in range(1, n + 1))
                whole = n in sizes
                ordered = all(m[i] >= m[i + 1] for i in range(q - 1))
                assert all_sizes == whole == ordered == hall(q, edges, m), (q, m)
                if n <= 6:
                    assert whole == perfect_matching(q, edges, m)
                total_weight += weight
                ordered_weight += weight * ordered
                distinct_weight += weight * (len(set(m)) == q)
                ordered_distinct_weight += weight * (ordered and len(set(m)) == q)
                cases += 1
            assert total_weight == q ** n
            assert factorial(q) * ordered_distinct_weight == distinct_weight
            if q == 2:
                tie = factorial(n) // factorial(n // 2) ** 2 if n % 2 == 0 else 0
                assert Q(ordered_weight, q ** n) == Q(1, 2) + Q(tie, 2 * q ** n)
            if n == 10:
                print(f'EXACT_FINITE_LAW: q={q}; n={n}; p_supported=1; Pr(N)={Q(ordered_weight,q**n)}; tie_free_order_fraction=1/{factorial(q)} when nonempty')
    print(f'ALL_SIZES: {cases} count vectors; q=2..5; n=1..10')


def zero_cases():
    # Lambda=0: even the directed q-cycle has no singleton cover.
    for q in range(2, 6):
        edges = frozenset((i, (i + 1) % q) for i in range(q))
        for i in range(q):
            assert not perfect_matching(q, edges, tuple(int(j == i) for j in range(q)))
    # Strongly connected looped triangle: P=conv(e1,(1,1,1)/3), dimension 1<2.
    edges = frozenset({(0, 0), (0, 1), (1, 2), (2, 0)})
    assert strongly_connected(3, edges)
    assert simple_cycles(3, edges) == {(0,), (0, 1, 2)}
    for n in range(1, 13):
        for m in compositions(n, 3):
            assert hall(3, edges, m) == (m[1] == m[2] and m[0] >= m[1])
    print('ZERO_CASES: no-loop singleton obstruction; dimension-deficient looped triangle forces N2=N3')


def main():
    source_fidelity()
    family_geometry()
    all_sizes_and_probability()
    zero_cases()
    print('SCOPE: exact finite checks and Gaussian chamber algebra; not a proof of the unbounded probabilistic limit')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
