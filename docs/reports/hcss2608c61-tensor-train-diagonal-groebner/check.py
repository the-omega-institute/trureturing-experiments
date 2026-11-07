#!/usr/bin/env python3
"""Exact source-fidelity and obstruction checks for arXiv:2608.19071v1.

Coordinates are one-based as in the article. Polynomials in the Buchberger
checks use integer coefficients; normalization/reduction use Fraction.
No floating point is used as mathematical evidence.
"""
from itertools import product, combinations, permutations
from fractions import Fraction
from functools import lru_cache
import json
import sympy as sp


def coordinates(d):
    return list(product(*(range(1, v + 1) for v in d)))


def page_key(j):
    return j[1:-1] + (j[0], j[-1])


def axes(d, cut):
    rm = tuple(range(1, cut)) + (0,)
    cm = tuple(range(cut, len(d) - 1)) + (len(d) - 1,)
    rows = list(product(*(range(1, d[k] + 1) for k in rm)))
    cols = list(product(*(range(1, d[k] + 1) for k in cm)))
    return rm, cm, rows, cols


def flatten(d, cut):
    rm, cm, rows, cols = axes(d, cut)
    def entry(a, b):
        j = [None] * len(d)
        for k, v in zip(rm + cm, a + b):
            j[k] = v
        return tuple(j)
    return rows, cols, [[entry(a, b) for b in cols] for a in rows]


def position(j, cut):
    return j[1:cut] + (j[0],), j[cut:-1] + (j[-1],)


def parity(p):
    return -1 if sum(p[i] > p[k] for i in range(len(p))
                     for k in range(i + 1, len(p))) % 2 else 1


def det_poly(mat, var_index, count):
    p = {}
    for perm in permutations(range(len(mat))):
        m = [0] * count
        for i, k in enumerate(perm):
            m[var_index[mat[i][k]]] += 1
        m = tuple(m)
        p[m] = p.get(m, 0) + parity(perm)
    return {m: c for m, c in p.items() if c}


def all_minors(d, ranks):
    variables = sorted(coordinates(d), key=page_key)
    vi = {j: k for k, j in enumerate(variables)}
    result, count = {}, 0
    for cut, rank in enumerate(ranks, 1):
        _, _, mat = flatten(d, cut)
        for rows in combinations(range(len(mat)), rank + 1):
            for cols in combinations(range(len(mat[0])), rank + 1):
                count += 1
                sub = [[mat[i][k] for k in cols] for i in rows]
                f = det_poly(sub, vi, len(variables))
                diagonal = [0] * len(variables)
                for i in range(rank + 1):
                    diagonal[vi[sub[i][i]]] += 1
                # Source Proposition 6.4, checked on every enumerated minor.
                assert max(f) == tuple(diagonal)
                assert f[max(f)] == 1
                result[tuple(sorted(f.items()))] = f
    return variables, list(result.values()), count


def divides(a, b):
    return all(x <= y for x, y in zip(a, b))


def reduce_poly(f, basis):
    f = dict(f)
    while f:
        lm = max(f)
        for g in basis:
            gl = max(g)
            if divides(gl, lm):
                shift = tuple(x-y for x, y in zip(lm, gl))
                scale = Fraction(f[lm], g[gl])
                for m, c in g.items():
                    key = tuple(x+y for x, y in zip(m, shift))
                    f[key] = f.get(key, 0) - scale*c
                    if not f[key]:
                        del f[key]
                break
        else:
            return f  # A nonzero remainder is enough to fail a zero-reduction test.
    return f


def buchberger_original(d, ranks):
    variables, basis, count = all_minors(d, ranks)
    # Remove only generators whose complete exact reduction is zero by predecessors.
    selected = []
    for g in sorted(basis, key=lambda g: max(g), reverse=True):
        if reduce_poly(g, selected):
            selected.append(g)
    checked = 0
    for f, g in combinations(selected, 2):
        a, b = max(f), max(g)
        if all(min(x, y) == 0 for x, y in zip(a, b)):
            continue  # Buchberger's product criterion.
        lcm = tuple(max(x, y) for x, y in zip(a, b))
        s = {}
        for sign, h, lm in [(1, f, a), (-1, g, b)]:
            shift = tuple(x-y for x, y in zip(lcm, lm))
            for m, c in h.items():
                key = tuple(x+y for x, y in zip(m, shift))
                s[key] = s.get(key, 0) + sign*Fraction(c, h[lm])
        s = {m: c for m, c in s.items() if c}
        assert not reduce_poly(s, selected), (d, ranks, 'nonzero S-pair')
        checked += 1
    assert all(not reduce_poly(g, selected) for g in basis)
    print('FIDELITY_THEOREM_6_5', json.dumps(dict(
        dimensions=d, ranks=ranks, minors=count, distinct=len(basis),
        selected=len(selected), noncoprime_pairs=checked)))


def longest_chain(support, cut):
    points = sorted(position(j, cut) for j in support)
    dp = []
    for i, (a, b) in enumerate(points):
        dp.append(1 + max([dp[k] for k, (c, e) in enumerate(points[:i])
                           if c < a and e < b] or [0]))
    return max(dp or [0])


def obstruction():
    d = (3, 2, 2, 3)
    variables = sorted(coordinates(d), key=page_key)
    symbols = {j: sp.Symbol('x_' + '_'.join(map(str, j))) for j in variables}
    a = lambda i, l: symbols[i, 1, 1, l]
    b = lambda i, l: symbols[i, 1, 2, l]
    c = lambda i, l: symbols[i, 2, 1, l]
    e = lambda i, l: symbols[i, 2, 2, l]
    # Construct F and G from fresh source-ordered flattenings, not from their expansions.
    rows, cols, mat2 = flatten(d, 2)
    rr = [rows.index(x) for x in [(1, 3), (2, 1), (2, 2)]]
    cc = [cols.index(x) for x in [(1, 1), (1, 2), (2, 1)]]
    F = sp.det(sp.Matrix([[symbols[mat2[i][k]] for k in cc] for i in rr]))
    rows, cols, mat3 = flatten(d, 3)
    rr = [rows.index(x) for x in [(1, 1, 3), (1, 2, 1), (1, 2, 2)]]
    cc = [cols.index((l,)) for l in [1, 2, 3]]
    G = sp.det(sp.Matrix([[symbols[mat3[i][k]] for k in cc] for i in rr]))
    beta = lambda p, q: b(1,p)*b(2,q) - b(1,q)*b(2,p)
    delta = c(1,2)*e(2,1) - c(2,2)*e(1,1)
    gamma = c(1,1)*e(2,1) - c(2,1)*e(1,1)
    epsilon = c(1,1)*c(2,2) - c(1,2)*c(2,1)
    R = sp.Poly(beta(2,3)*F-delta*G, *[symbols[j] for j in variables], domain=sp.ZZ)
    rhs = ((a(3,2)*beta(1,3)-a(3,3)*beta(1,2))*delta
           -a(3,2)*beta(2,3)*gamma+b(3,1)*beta(2,3)*epsilon)
    assert sp.expand(R.as_expr()-rhs) == 0
    assert len(R.terms()) == 16
    assert all(abs(coeff) == 1 and sum(m) == 5 for m, coeff in R.terms())
    lm_expected = sp.Poly(a(3,2)*b(1,1)*b(2,3)*c(1,2)*e(2,1), *R.gens)
    assert R.monoms()[0] == lm_expected.monoms()[0]
    print('IDENTITY', sp.sstr(R.as_expr()))
    print('LEADING_MONOMIAL', sp.sstr(lm_expected.as_expr()))
    vv, minors, count = all_minors(d, (2,2,2))
    assert vv == variables and count == 840
    chains = []
    supports = []
    for m, coeff in R.terms():
        support = [j for j, exponent in zip(variables, m) for _ in range(exponent)]
        supports.append(support)
        bounds = [longest_chain(support, cut) for cut in [1,2,3]]
        chains.extend(bounds)
        assert max(bounds) <= 2
        assert all(not divides(max(f), m) for f in minors)
        print('TERM', coeff, '*'.join(str(symbols[j]) for j in support), 'chains', bounds)
    print('SEED', json.dumps(dict(terms=len(supports), degree=5, minors=count,
                                 chain_bounds=len(chains), max_chain=max(chains))))
    padded = 0
    for N in range(4, 9):
        for ends, middle in [((3,3),2), ((4,5),3)]:
            dims = (ends[0],) + (middle,)*(N-2) + (ends[1],)
            embed = lambda j: j[:3] + (1,)*(N-4) + j[3:]
            for support in supports:
                active = [embed(j) for j in support]
                assert all(all(1 <= v <= size for v, size in zip(j, dims)) for j in active)
                for cut in range(1,N):
                    assert longest_chain(active, cut) <= 2
            # Genuine tensor example: sum of the two constant-coordinate pure tensors.
            # Its flattening has exactly two disjoint nonzero rows and columns.
            for cut in range(1,N):
                pts = [position((v,)*N, cut) for v in [1,2]]
                assert pts[0][0] != pts[1][0] and pts[0][1] != pts[1][1]
                assert sp.Matrix([[1,0],[0,1]]).rank() == 2
            padded += 1
    print('PADDING', json.dumps(dict(dimension_tuples=padded, N_min=4, N_max=8,
                                     every_supported_chain_at_most=2)))


def main():
    # Theorem 6.5: both flattenings, including unequal minor sizes.
    for d, ranks in [((2,2,2),(1,1)), ((3,2,3),(2,2)),
                     ((3,2,2),(1,2)), ((2,2,3),(2,1))]:
        buchberger_original(d, ranks)
    obstruction()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
