#!/usr/bin/env python3
"""Exact preregistration checks for CKS arXiv:2610.37841v1, question (4)(a).

No Lean proof is produced. The finite checks support the general argument in
README.md; they do not certify the arbitrary-arity or infinite-family claims.
Python standard library only; all probability and interpolation data are rational.
"""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations, product
from math import gcd


S = ((0, 0, 0), (0, 2, 0), (0, 0, 0))  # ordered 0,a,b
Q = tuple(map(F, ('1/3', '5/9', '17/27', '2/3', '7/9',
                 '65/81', '23/27', '8/9', '1')))


def probability(left, right):
    assert len(left) == len(right) > 0
    return F(sum(x == y for x, y in zip(left, right)), len(left))


def source_fidelity():
    """Example 4.1, formula (star), respecting its positive-side-arity guards."""
    checked = 0
    for a, b, c in product(range(4), repeat=3):
        if a + c == 0 or b + c == 0:
            continue
        assignments = tuple(product(range(2), repeat=a+b+c))
        left = tuple(int(all(v[:c+a])) for v in assignments)
        right = tuple(int(all(v[:c]) and all(v[c+a:])) for v in assignments)
        expected = 1 - F(2**a + 2**b - 2, 2**(a+b+c))
        assert probability(left, right) == expected, (a, b, c)
        checked += 1
    assert 1 - F(2**1 + 2**1 - 2, 2**3) == F(3, 4)
    print(f'FIDELITY: Example 4.1 formula (star), {checked} guarded triples in [0,3]^3; xy=xz gives 3/4')


def normal_forms(k):
    assignments = tuple(product(range(3), repeat=k))
    forms = {'zero': (0,) * len(assignments)}
    forms.update({f'x{i}': tuple(v[i] for v in assignments) for i in range(k)})
    forms.update({f'x{i}^2': tuple(S[v[i]][v[i]] for v in assignments) for i in range(k)})
    forms.update({f'x{i}x{j}': tuple(S[v[i]][v[j]] for v in assignments)
                  for i, j in combinations(range(k), 2)})
    return assignments, forms


def term_closure(k):
    """Generate the actual clone component from projections and binary multiplication."""
    assignments = tuple(product(range(3), repeat=k))
    tables = {tuple(v[i] for v in assignments) for i in range(k)}
    while True:
        old = tuple(tables)
        tables.update(tuple(S[x][y] for x, y in zip(f, g)) for f in old for g in old)
        if len(tables) == len(old):
            return tables


def small_semigroup():
    for x, y, z in product(range(3), repeat=3):
        assert S[S[x][y]][z] == S[x][S[y][z]] == 0
    assert all(S[x][y] == S[y][x] for x, y in product(range(3), repeat=2))
    for k in range(1, 7):
        _, forms = normal_forms(k)
        generated = term_closure(k)
        assert generated == set(forms.values())
        observed = {probability(f, g) for f in generated for g in generated}
        assert observed <= set(Q)
        if k >= 4:
            assert observed == set(Q)
        print(f'S: arity {k}, generated term functions {len(generated)}, spectrum values '+','.join(map(str, sorted(observed))))
    _, forms = normal_forms(4)
    cases = [('zero', 'x0^2', F(2, 3)), ('zero', 'x0x1', F(8, 9)),
             ('x0^2', 'x1^2', F(5, 9)), ('x0^2', 'x0x1', F(7, 9)),
             ('x0^2', 'x1x2', F(17, 27)), ('x0x1', 'x0x2', F(23, 27)),
             ('x0x1', 'x2x3', F(65, 81)), ('x0', 'x0^2', F(1, 3)),
             ('x0', 'x1x2', F(1, 3)), ('x0', 'x1', F(1, 3)),
             ('x0x1', 'x0x1', F(1))]
    for left, right, expected in cases:
        assert probability(forms[left], forms[right]) == expected
    assert S[0][0] == S[0][1] == S[1][0] == 0 and S[1][1] == 2
    print('S: all eleven comparison types checked; term-condition failure t(x,y)=xy at a=0,b=a,u=0,v=a')


def cyclic_checks():
    for e in range(1, 9):
        n = 3**e
        coefficients = range(n) if e <= 5 else sorted({0, 1, 2, n-1} | {3**i for i in range(e)})
        tested = 0
        for c in coefficients:
            fibers = Counter(c*x % n for x in range(n))
            g = gcd(n, c)
            assert len(fibers) == n // g and set(fibers.values()) == {g}
            assert F(fibers[0], n) in {F(1, 3**j) for j in range(e+1)}
            tested += 1
        print(f'C: modulus 3^{e}, {tested} unary coefficient residues, exact uniform fibers')
    checked = 0
    for e in range(1, 4):
        n = 3**e
        for c, d in product(range(n), repeat=2):
            fibers = Counter((c*x+d*y) % n for x, y in product(range(n), repeat=2))
            g = gcd(n, gcd(c, d))
            assert len(fibers) == n // g and set(fibers.values()) == {n*g}
            assert F(fibers[0], n*n) == F(g, n)
            checked += 1
    print(f'C: {checked} complete binary coefficient vectors over moduli 3,9,27')
    n = 3**8
    attained = []
    for j in range(9):
        coefficient = 3**(8-j)+1
        p = F(sum((coefficient*x) % n == x for x in range(n)), n)
        assert p == F(1, 3**j)
        attained.append(p)
    print('C=Z/6561Z: all nine native equations (3^(8-j)+1)x=x attain '+','.join(map(str, attained)))
    checked = 0
    for n in (3, 9):
        for c, d in product(range(n), repeat=2):
            for a, b, u, v in product(range(n), repeat=4):
                if (c*a+d*u) % n == (c*a+d*v) % n:
                    assert (c*b+d*u) % n == (c*b+d*v) % n
                checked += 1
    print(f'C: {checked} one-parameter linear term-condition instances over Z/3Z and Z/9Z')


def finite_power_checks():
    checks = 0
    for m, k in ((2, 1), (2, 2), (2, 3), (3, 1), (3, 2)):
        assignments, forms = normal_forms(k)
        table_values = tuple(forms.values())
        powered_assignments = tuple(product(range(len(assignments)), repeat=m))
        for f, g in product(table_values, repeat=2):
            powered_count = sum(all(f[i] == g[i] for i in v) for v in powered_assignments)
            assert F(powered_count, len(powered_assignments)) == probability(f, g)**m
            checks += 1
    print(f'POWER LAW: {checks} native normal-form equations, direct products of 2 or 3 S coordinates')


def interpolate(x, xs, ys):
    assert xs[0] <= x <= xs[-1]
    for i in range(len(xs)-1):
        if x <= xs[i+1]:
            return ys[i] + (x-xs[i]) * (ys[i+1]-ys[i]) / (xs[i+1]-xs[i])
    raise AssertionError('uncovered endpoint')


def family_checks():
    assert all(x < y for x, y in zip(Q, Q[1:]))
    for d in range(1, 7):
        left = (F(0),) + tuple(q**(8*d) for q in Q)
        right = (F(0),) + tuple(F(1, 3**(j*d)) for j in reversed(range(9)))
        assert len(set(left)) == len(set(right)) == 10
        assert left[1] == right[1] == F(1, 3**(8*d))
        assert left[-1] == right[-1] == 1
        assert 3**(8*d) == (3**8)**d
        for i in range(9):
            slope = (right[i+1]-right[i]) / (left[i+1]-left[i])
            assert slope > 0
            midpoint = (left[i]+left[i+1]) / 2
            assert interpolate(interpolate(midpoint, left, right), right, left) == midpoint
        assert all(interpolate(x, left, right) == y for x, y in zip(left, right))
        zeros = (0,) * (8*d)
        a = (1,) + zeros[1:]
        mul = lambda x, y: tuple(S[u][v] for u, v in zip(x, y))
        assert mul(zeros, zeros) == mul(zeros, a) == mul(a, zeros) == zeros
        assert mul(a, a) == (2,) + zeros[1:]
        print(f'FAMILY: d={d}, common cardinality {3**(8*d)}, nine spectral knots and nine positive slopes; inverse and nonabelian witness checked')


def main():
    source_fidelity()
    small_semigroup()
    cyclic_checks()
    finite_power_checks()
    family_checks()
    print('SCOPE: finite checks above; the arbitrary-arity normal-form induction, general linear cancellation, power law and all-d homeomorphism argument remain proposed proofs')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
