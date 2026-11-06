#!/usr/bin/env python3
"""Exact Rao-cokernel checks for arXiv:2605.19434v1, Question 6.3.

For distinct u_i in P^1(Q), parametrize the lines by
(z0,z1,z2,z3)=(u_i0*v0,u_i0*v1,u_i1*v0,u_i1*v1).
The quadric is z0*z3-z1*z2. Proposition 4.2 identifies M_t with the
cokernel of R_t -> direct_sum_i Q[v0,v1]_t, compatibly with multiplication.
All ranks below are over Q, computed with SymPy exact integer matrices.
No finite collection of these checks proves the unbounded question.
"""
from functools import lru_cache
from math import comb
import sys
import sympy as sp


def compositions(n, parts):
    if parts == 1:
        yield (n,)
    else:
        for k in range(n + 1):
            for tail in compositions(n - k, parts - 1):
                yield (k,) + tail


@lru_cache(None)
def restriction(points, t):
    """Columns are restrictions of all ambient degree-t monomials."""
    if t < 0:
        return sp.zeros(0, 0)
    exps = list(compositions(t, 4))
    out = sp.zeros(len(points) * (t + 1), len(exps))
    for col, e in enumerate(exps):
        for i, (a, b) in enumerate(points):
            out[i * (t + 1) + e[1] + e[3], col] = (
                a ** (e[0] + e[1]) * b ** (e[2] + e[3]))
    return out


@lru_cache(None)
def restriction_basis(points, t):
    mat = restriction(points, t)
    basis = mat.columnspace()
    return sp.Matrix.hstack(*basis) if basis else sp.zeros(mat.rows, 0)


def dim_rao(points, t):
    if t < 0:
        return 0
    return len(points) * (t + 1) - restriction_basis(points, t).cols


def multiply(points, t, m, c):
    """Actual multiplication by (sum c_j*z_j)^m on the line sections."""
    if t < 0:
        return sp.zeros(len(points) * max(t + m + 1, 0), 0)
    out = sp.zeros(len(points) * (t + m + 1), len(points) * (t + 1))
    for i, (a, b) in enumerate(points):
        A, B = c[0] * a + c[2] * b, c[1] * a + c[3] * b
        for j in range(t + 1):
            for k in range(m + 1):
                out[i * (t + m + 1) + j + k, i * (t + 1) + j] = (
                    comb(m, k) * A ** (m - k) * B ** k)
    return out


@lru_cache(None)
def rank_rao(points, t, m, c):
    if t < 0:
        return 0
    target = restriction_basis(points, t + m)
    mult = multiply(points, t, m, c)
    # Well-definedness on the actual quotient, checked independently.
    source = restriction_basis(points, t)
    assert target.row_join(mult * source).rank() == target.cols
    return target.row_join(mult).rank() - target.cols


def differential(d, t, m):
    """Independent proposed model E^m, with ordinary monomial bases."""
    if not 0 <= t <= t + m <= d:
        return sp.zeros(0, 0)
    src = [(a, d - t - a, b, t - b)
           for a in range(d - t + 1) for b in range(t + 1)]
    dst = [(a, d - t - m - a, b, t + m - b)
           for a in range(d - t - m + 1) for b in range(t + m + 1)]
    ix = {e: i for i, e in enumerate(dst)}
    out = sp.zeros(len(dst), len(src))
    for j, (a, b, c, e) in enumerate(src):
        for k in range(m + 1):
            ell = m - k
            if k <= a and ell <= b:
                coeff = (comb(m, k) * sp.factorial(a) // sp.factorial(a-k)
                         * sp.factorial(b) // sp.factorial(b-ell))
                out[ix[a-k, b-ell, c+k, e+ell], j] += coeff
    return out


def check(points):
    r = len(points)
    assert all(a*d-b*c != 0 for i, (a,b) in enumerate(points)
               for c,d in points[i+1:])
    # Source fidelity: Lemma 4.1 (symmetric Hilbert function), Proposition
    # 4.2 (support 0..r-2; dim M_0=r-1), and Proposition 4.4 (WLP).
    support = [dim_rao(points, t) for t in range(r + 1)]
    expected = [max((t + 1) * (r - 1 - t), 0) for t in range(r + 1)]
    assert support == expected
    assert support[:r-1] == support[:r-1][::-1]
    assert support[0] == r - 1
    good = [(1, 0, 0, 1), (2, 3, 5, 7)]
    slp_count = 0
    for c in good:
        assert c[0]*c[3]-c[1]*c[2] != 0
        for t in range(r - 1):
            for m in range(1, r - t + 1):
                got = rank_rao(points, t, m, c)
                want = min(dim_rao(points, t), dim_rao(points, t + m))
                assert got == want, (r, points, c, t, m, got, want)
                slp_count += 1
    # Extract Jordan block multiplicities from successive global ranks.
    d = r - 2
    ranks = [sum(support)] + [
        sum(rank_rao(points, t, m, good[0]) for t in range(max(d - m + 1, 0)))
        for m in range(1, d + 3)]
    blocks = {k: ranks[k-1]-2*ranks[k]+ranks[k+1]
              for k in range(1, d+2) if ranks[k-1]-2*ranks[k]+ranks[k+1]}
    predicted = {d-2*j+1: d-2*j+1 for j in range(d//2+1)}
    assert blocks == predicted
    for t in range(d + 1):
        for m in range(1, d - t + 1):
            assert differential(d, t, m).rank() == rank_rao(points, t, m, good[0])
    if r >= 3:
        c = (1, 2, 2, 4)  # rank-one coefficient matrix: a tangent plane
        got = rank_rao(points, 0, 1, c)
        want = min(dim_rao(points, 0), dim_rao(points, 1))
        assert got == d and got < want
        assert rank_rao(points, 0, 1, (0, 0, 0, 0)) == 0
        print(f'NEGATIVE_CONTROL r={r} tangent_rank={got} maximal_rank={want}')
    print(f'r={r} points={points} hilbert={support} '
          f'power_maps={slp_count} jordan_length_multiplicity={blocks}')
    return slp_count


def main():
    print('Python', sys.version.split()[0], 'SymPy', sp.__version__, 'field=Q exact')
    total = 0
    for r in range(2, 7):
        for points in [tuple((1, a) for a in range(r)),
                       ((0, 1),) + tuple((1, a) for a in range(r-1))]:
            total += check(points)
    one_line = ((1, 0),)
    assert all(dim_rao(one_line, t) == 0 for t in range(5))
    print('FIDELITY_OK: source Hilbert symmetry, initial dimension, support, WLP')
    print(f'SLP_OK: {total} exact induced power maps on 10 actual line configurations')
    print('MODEL_OK: differential ranks and predicted Jordan blocks agree')
    print('BOUNDARY_OK: r=1 tested; tangent and zero-form negative controls detected')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
