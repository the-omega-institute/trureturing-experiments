#!/usr/bin/env python3
"""Exact monomial and Taylor-resolution checks for arXiv:2609.37219, Q4.10.

No regularity bound is assumed in taylor_betti: it tensors the Taylor free
resolution of an ideal with the residue field, then computes homology by exact
row reduction. Characteristic zero uses Fraction; the other fields are F_2,F_3.
All exponent windows and parameter ranges are printed. These finite checks do
not prove the quantified question or the homological bridges needed by Lean.
"""
from collections import defaultdict
from fractions import Fraction
from itertools import combinations, permutations, product


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def divides(a, b):
    return all(x <= y for x, y in zip(a, b))


def edges(m):
    return [(u, v) for u in range(m*m) for v in range(u+1, m*m)
            if u//m == v//m or u % m == v % m]


def generators(m):
    return {tuple(int(j != s[i]) for i in range(m) for j in range(m))
            for s in permutations(range(m))}


def ordinary(m, q):
    g = generators(m)
    current = {(0,) * (m*m)}
    for _ in range(q):
        current = {add(a, b) for a in current for b in g}
    return current


def symbolic_minimal(m, q):
    ee = edges(m)
    answer = set()
    for a in product(range(q+1), repeat=m*m):
        if not all(a[u] + a[v] >= q for u, v in ee):
            continue
        # Decreasing a positive coordinate by 1 must violate an incident edge.
        if all(x == 0 or any(a[u]+a[v] == q for u, v in ee if k in (u, v))
               for k, x in enumerate(a)):
            answer.add(a)
    return answer


def rank(matrix, characteristic):
    if not matrix or not matrix[0]:
        return 0
    if characteristic:
        a = [[x % characteristic for x in row] for row in matrix]
    else:
        a = [[Fraction(x) for x in row] for row in matrix]
    r = 0
    for col in range(len(a[0])):
        pivot = next((i for i in range(r, len(a)) if a[i][col]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        inv = pow(a[r][col], -1, characteristic) if characteristic else 1 / a[r][col]
        a[r] = [v * inv % characteristic if characteristic else v * inv for v in a[r]]
        for i in range(len(a)):
            if i == r or not a[i][col]:
                continue
            c = a[i][col]
            a[i] = [(v-c*w) % characteristic if characteristic else v-c*w
                    for v, w in zip(a[i], a[r])]
        r += 1
        if r == len(a):
            break
    return r


def taylor_betti(gens, characteristic):
    gg = sorted(gens)
    groups = defaultdict(lambda: defaultdict(list))
    lcms = {}
    for size in range(1, len(gg)+1):
        for subset in combinations(range(len(gg)), size):
            lcm = tuple(max(gg[k][v] for k in subset) for v in range(len(gg[0])))
            lcms[subset] = lcm
            groups[lcm][size-1].append(subset)
    betti = defaultdict(int)
    for lcm, chain in groups.items():
        ranks = {}
        for i, source in chain.items():
            target = chain.get(i-1, [])
            rows = {s: k for k, s in enumerate(target)}
            mat = [[0]*len(source) for _ in target]
            for col, subset in enumerate(source):
                for position in range(len(subset)):
                    deleted = subset[:position] + subset[position+1:]
                    if deleted in rows:
                        mat[rows[deleted]][col] = (-1) ** position
            ranks[i] = rank(mat, characteristic)
        for i, basis in chain.items():
            dimension = len(basis) - ranks.get(i, 0) - ranks.get(i+1, 0)
            assert dimension >= 0
            if dimension:
                betti[(i, sum(lcm))] += dimension
    return dict(sorted(betti.items()))


def check_regularities():
    for m, qrange in [(1, range(3, 7)), (2, range(1, 7))]:
        for q in qrange:
            og = ordinary(m, q)
            sg = symbolic_minimal(m, q)
            assert sg == og
            by_field = []
            for p in (0, 2, 3):
                b = taylor_betti(sg, p)
                expected = {(0, 0): 1} if m == 1 else {(0, 2*q): q+1, (1, 2*q+2): q}
                assert b == expected, (m, q, p, b)
                reg = max(j-i for (i, j), multiplicity in b.items() if multiplicity)
                assert reg == (0 if m == 1 else 2*q+1)
                by_field.append((p, reg))
            print('EXACT_BETTI m=%d q=%d beta=%s fields_and_reg=%s' % (m, q, b, by_field))
    print('FIDELITY_OK: Theorem 4.7 at m=n=2,q=2..6; Corollary 4.9 at m=n=2,q=1..6')
    print('SMALL_SETTLEMENT_OK: Q4.10 at m=1,2,q=3..6, over Q,F_2,F_3')


def check_isolated_lcm():
    for m in range(2, 5):
        U = tuple(int(i != j) for i in range(m) for j in range(m))
        s = list(range(m)); s[0], s[1] = s[1], s[0]
        V = tuple(int(j != s[i]) for i in range(m) for j in range(m))
        for q in range(1, 5):
            gg = ordinary(m, q)
            P = tuple(q*x for x in U)
            Q = add(tuple((q-1)*x for x in U), V)
            L = tuple(max(x, y) for x, y in zip(P, Q))
            below = {a for a in gg if divides(a, L)}
            assert below == {P, Q}
            assert sum(L) == q*m*(m-1)+2
            for k, x in enumerate(L):
                if x:
                    lower = list(L); lower[k] -= 1
                    assert sum(divides(a, lower) for a in below) <= 1
            print('ISOLATED_LCM m=%d q=%d distinct_generators=%d below_L=2 degree_L=%d' %
                  (m, q, len(gg), sum(L)))


def check_colons():
    total = 0
    for m in (2, 3):
        ee = edges(m); universe = set(range(m*m))
        for q in range(2, 5):
            checked = 0; independent_sets = 0
            for mask in range(1 << (m*m)):
                A = {u for u in universe if (mask >> u) & 1}
                if any(u in A and v in A for u, v in ee):
                    # For A containing an edge every symbolic monomial lies in (A).
                    continue
                independent_sets += 1
                B = sorted(universe - A)
                N = {v if u in A else u for u, v in ee if (u in A) != (v in A)}
                residual = universe - A - N
                a = len(A); r = m-a
                assert len(N) == a*(2*m-a-1)
                assert len(B) == m*(m-1)+r
                assert len(residual) == r*r
                for values in product(range(q), repeat=len(B)):
                    e = dict(zip(B, values))
                    shifted = {u: 0 if u in A else e[u]+1 for u in universe}
                    lhs = all(shifted[u]+shifted[v] >= q for u, v in ee)
                    rhs = (all(e[v] >= q-1 for v in N) and
                           all(e[u]+e[v] >= q-2 for u, v in ee
                               if u in residual and v in residual))
                    assert lhs == rhs, (m, q, A, e)
                    checked += 1
                # Monomials containing a variable of A belong to both ideals.
            total += checked
            print('COLON_IDENTITY m=%d q=%d independent_sets=%d exponent_window=0..%d cases=%d' %
                  (m, q, independent_sets, q-1, checked))
    print('COLON_CASES_TOTAL', total)
    # Exact integer form of the induction's polynomial bounds, checked on a window.
    for m in range(2, 31):
        d = m*(m-1)
        for a in range(m+1):
            r = m-a; dr = r*(r-1)
            for q in range(2, 31):
                R = 0 if r <= 1 or q == 2 else (r*r-1 if q == 3 else (q-2)*dr+1)
                contribution = (q-1)*(d-dr)+R+d+r
                assert contribution == q*d+(R-(q-1)*dr+r)
                assert contribution <= q*d+1
    print('INDUCTION_ARITHMETIC_OK: 2<=m,q<=30,0<=a<=m; exact integers')


def check_source_caution():
    # The source's separate m=1 rectangular proof asserts a false equality.
    triangle = [(0, 1), (0, 2), (1, 2)]
    xyz = (1, 1, 1)
    assert all(xyz[u]+xyz[v] >= 2 for u, v in triangle)
    minimal_cover = {(1, 1, 0), (1, 0, 1), (0, 1, 1)}
    square = {add(a, b) for a in minimal_cover for b in minimal_cover}
    assert not any(divides(g, xyz) for g in square)
    print('SOURCE_CAUTION_CONFIRMED: xyz in J(K_3)^(2) but not J(K_3)^2; unused by square induction')


def main():
    check_regularities()
    check_isolated_lcm()
    check_colons()
    check_source_caution()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
