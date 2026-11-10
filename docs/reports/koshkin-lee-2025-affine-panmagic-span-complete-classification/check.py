#!/usr/bin/env python3
"""Exact finite checks for the Koshkin--Lee affine-span question.

Self-contained Python standard library; integer row operations compute Q-ranks.
Finite checks support, but do not prove, the uniform argument in the draft.
"""

from itertools import permutations
from math import gcd, prod


def factor(n):
    ans = []
    p = 2
    while p * p <= n:
        e = 0
        while n % p == 0:
            n //= p
            e += 1
        if e:
            ans.append((p, e))
        p += 1
    if n > 1:
        ans.append((n, 1))
    return ans


def slopes(n):
    return [a for a in range(n) if all(gcd(a + d, n) == 1 for d in (-1, 0, 1))]


def rank_formula(n):
    return prod(1 + (p - 1) * (p - 3) * sum(p ** (2 * j) for j in range(e))
                for p, e in factor(n))


def fourier_support(n):
    return {((-a * s) % n, s) for a in slopes(n) for s in range(n)}


def panmagic_support(n):
    return {(0, 0)} | {(r, s) for r in range(n) for s in range(n)
                        if r and s and (r + s) % n and (r - s) % n}


def rational_rank(rows):
    """Exact elimination: nonzero integer scaling, row subtraction, gcd reduction."""
    basis = {}
    for row in rows:
        v = {i: x for i, x in enumerate(row) if x}
        while v:
            p = min(v)
            if p not in basis:
                g = gcd(*v.values())
                sign = 1 if v[p] > 0 else -1
                basis[p] = {i: sign * x // g for i, x in v.items()}
                break
            b = basis[p]
            vp, bp = v[p], b[p]
            v = {i: bp * v.get(i, 0) - vp * b.get(i, 0) for i in v.keys() | b.keys()}
            v = {i: x for i, x in v.items() if x}
            if v:
                g = gcd(*v.values())
                v = {i: x // g for i, x in v.items()}
    return len(basis)


def permutation_rows(n):
    # Transposed source convention, with input x first: y = a*x+b.
    for a in slopes(n):
        for b in range(n):
            yield [int(y == (a * x + b) % n) for x in range(n) for y in range(n)]


def line_sums(m):
    n = len(m)
    return ([sum(m[x][y] for y in range(n)) for x in range(n)]
            + [sum(m[x][y] for x in range(n)) for y in range(n)]
            + [sum(m[x][(k - x) % n] for x in range(n)) for k in range(n)]
            + [sum(m[x][(x - k) % n] for x in range(n)) for k in range(n)])


def target_dimension(n):
    lines = []
    for family in range(4):
        for k in range(n):
            lines.append([int((x == k, y == k, (x + y) % n == k,
                               (x - y) % n == k)[family])
                          for x in range(n) for y in range(n)])
    reference = lines[0]
    return n * n - rational_rank(([x - y for x, y in zip(line, reference)]
                                  for line in lines[1:]))


def mobius(n):
    fs = factor(n)
    return 0 if any(e > 1 for _, e in fs) else (-1) ** len(fs)


def ramanujan(n, k):
    return sum(d * mobius(n // d) for d in range(1, n + 1)
               if n % d == 0 and k % d == 0)


def check_witness(n):
    p = factor(n)[0][0]
    phi = sum(gcd(t, n) == 1 for t in range(n))
    table = [ramanujan(n, k) for k in range(n)]
    m = [[table[(p * x + y) % n] for y in range(n)] for x in range(n)]
    assert m[0][0] == phi > 0
    assert set(line_sums(m)) == {0}
    assert sum(x * x for row in m for x in row) == n * n * phi
    for a in slopes(n):
        for b in range(n):
            assert sum(m[x][(a * x + b) % n] for x in range(n)) == 0
    # epsilon = 1/(2*n*phi): scaled positive matrix has entries 2*phi+R.
    positive = [[2 * phi + x for x in row] for row in m]
    assert min(min(row) for row in positive) > 0
    assert set(line_sums(positive)) == {2 * n * phi}
    return p, phi


def main():
    # Fidelity: source Section 4's ten order-five generators have exactly one relation.
    n = 5
    pan_perms = [q for q in permutations(range(n))
                 if len({(q[x] - x) % n for x in range(n)}) == n
                 and len({(q[x] + x) % n for x in range(n)}) == n]
    affine_perms = {tuple((a * x + b) % n for x in range(n))
                    for a in slopes(n) for b in range(n)}
    assert set(pan_perms) == affine_perms and len(affine_perms) == 10
    assert rational_rank(permutation_rows(5)) == target_dimension(5) == 9
    assert all(sum(int(y == (a * x + b) % 5) for b in range(5)) == 1
               for a in (2, 3) for x in range(5) for y in range(5))
    print('FIDELITY n=5: all 10 panmagic permutations affine; exact Q-rank=9; target dimension=9; sum_b P_2,b=sum_b P_3,b=J')
    # Corollary 1 (journal Corollary 5.2), checked from literal diagonal counts.
    for n in range(2, 61):
        allowed = slopes(n)
        for a in range(n):
            if gcd(a, n) == 1:
                for b in range(n):
                    actual = (len({((a - 1) * x + b) % n for x in range(n)}) == n
                              and len({((a + 1) * x + b) % n for x in range(n)}) == n)
                    assert actual == (a in allowed)
    print('FIDELITY affine panmagic iff a,a-1,a+1 units: every affine map at 2<=n<=60 checked')
    orders = [n for n in range(5, 151) if gcd(n, 6) == 1]
    composites = []
    print('n support_cardinality product_formula panmagic_dimension prime equality')
    for n in orders:
        support = fourier_support(n)
        target = panmagic_support(n)
        prime = factor(n) == [(n, 1)]
        assert support <= target
        assert len(support) == rank_formula(n)
        assert len(target) == (n - 2) ** 2
        assert (support == target) == prime
        # CRT support equality, tested as actual sets, including repeated prime powers.
        local_moduli = [p ** e for p, e in factor(n)]
        local_sets = [fourier_support(q) for q in local_moduli]
        crt_image = {tuple((r % q, s % q) for q in local_moduli) for r, s in support}
        assert len(crt_image) == prod(len(local) for local in local_sets)
        assert all(pair in local for point in crt_image
                   for local, pair in zip(local_sets, point))
        print(n, len(support), rank_formula(n), len(target), prime, support == target)
        if not prime:
            composites.append(n)
    print(f'SUPPORT_AND_CRT {len(orders)} Polya orders checked, including {len(composites)} composites')
    for n in (5, 7, 11, 13, 25, 35):
        rank = rational_rank(permutation_rows(n))
        dim = target_dimension(n)
        assert rank == rank_formula(n) and dim == (n - 2) ** 2
        print(f'EXACT_Q_RANK n={n} affine={rank} panmagic={dim}')
    for n in composites:
        p, phi = check_witness(n)
        print(f'INTEGER_WITNESS n={n} p={p} phi={phi} four_zero_line_families=True all_affine_inner_products_zero=True norm_squared={n*n*phi} positive_perturbation=True')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
