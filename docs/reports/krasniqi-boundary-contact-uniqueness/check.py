#!/usr/bin/env python3
"""Exact preregistration checks for Krasniqi Conjecture 23.

Fidelity: reconstruct source equations (3.17)-(3.19) as formal rational
polynomials, independently via the exponential generating product.
Expected-settlement check: verify every closed box of the repository's
uniform filtered-density certificate and its unbounded-tail premise.
This program does not evaluate the source improper integral or prove the
conjecture. The analytic bridges are stated in README.md.

All evidence uses integers and fractions.Fraction; no floating point.
Partition data: KrasniqiBoundaryContact.md, Appendix A, main-repository
revision 458aaeee2237c03eff58e0bb81b72d1d7a7431a1.
"""
from collections import Counter
from fractions import Fraction as Q
from functools import lru_cache
from math import comb, factorial

# Sparse polynomials in a,c,t with rational coefficients.
ZERO = {}
ONE = {(0, 0, 0): Q(1)}
A = {(1, 0, 0): Q(1)}
C = {(0, 1, 0): Q(1)}


def add(*polys):
    ans = {}
    for poly in polys:
        for e, v in poly.items():
            ans[e] = ans.get(e, Q(0)) + v
    return {e: v for e, v in ans.items() if v}


def scale(poly, v):
    return {e: w * v for e, w in poly.items() if w * v}


def mul(p, q):
    ans = {}
    for e, v in p.items():
        for f, w in q.items():
            g = tuple(x + y for x, y in zip(e, f))
            ans[g] = ans.get(g, Q(0)) + v * w
    return {e: v for e, v in ans.items() if v}


def power(p, n):
    ans = ONE
    for _ in range(n):
        ans = mul(ans, p)
    return ans


def exact_coefficients(n):
    qs = [ONE]
    for j in range(1, n + 1):
        qs.append(scale(add(*(mul(add(scale(A, Q(k, k + 1)), C), qs[j-k])
                              for k in range(1, j + 1))), Q(1, j)))
    # Independent truncated product of exp(ell_k z^k).
    product = [ONE] + [{} for _ in range(n)]
    for k in range(1, n + 1):
        ell = add(scale(A, Q(1, k + 1)), scale(C, Q(1, k)))
        old = product
        product = [{} for _ in range(n + 1)]
        for j in range(n + 1):
            for m in range(j // k + 1):
                product[j] = add(product[j], scale(mul(old[j-k*m], power(ell, m)),
                                                  Q(1, factorial(m))))
    assert product == qs
    return qs


def canonical(poly):
    return tuple(sorted(poly.items()))


@lru_cache(maxsize=None)
def matrix(n, lo, hi):
    """Power monomials on [lo,hi] to Bernstein basis of degree n."""
    return tuple(tuple(sum((Q(comb(j, u)) * lo ** (j-u) * (hi-lo) ** u *
                             Q(comb(i, u), comb(n, u))
                             for u in range(min(i, j) + 1)), Q(0))
                       for j in range(n + 1)) for i in range(n + 1))


@lru_cache(maxsize=None)
def transform(poly, axis, n, lo, hi):
    mat = matrix(n, lo, hi)
    ans = {}
    for e, v in poly:
        for i in range(n + 1):
            w = mat[i][e[axis]]
            if w:
                f = e[:axis] + (i,) + e[axis+1:]
                ans[f] = ans.get(f, Q(0)) + v * w
    return canonical({e: v for e, v in ans.items() if v})


def bernstein_values(poly, box):
    degrees = tuple(max(e[i] for e in poly) for i in range(3))
    value = canonical(poly)
    for i in range(3):
        value = transform(value, i, degrees[i], *box[i])
    values = dict(value)
    return [values.get((i, j, k), Q(0))
            for i in range(degrees[0] + 1)
            for j in range(degrees[1] + 1)
            for k in range(degrees[2] + 1)]


# Complete closed-face binary partition. 0=a, 1=c, 2=t.
LEAVES = """0L1L2L0L k
0L1L2L0R1L k
0L1L2L0R1R2L0L1L k
0L1L2L0R1R2L0L1R2L0L k
0L1L2L0R1R2L0L1R2L0R1L k
0L1L2L0R1R2L0L1R2L0R1R +
0L1L2L0R1R2L0L1R2R0L k
0L1L2L0R1R2L0L1R2R0R +
0L1L2L0R1R2L0R1L2L0L k
0L1L2L0R1R2L0R1L2L0R +
0L1L2L0R1R2L0R1L2R0L k
0L1L2L0R1R2L0R1L2R0R +
0L1L2L0R1R2L0R1R +
0L1L2L0R1R2R +
0L1L2R +
0L1R2L0L1L2L0L k
0L1R2L0L1L2L0R1L k
0L1R2L0L1L2L0R1R2L0L1L k
0L1R2L0L1L2L0R1R2L0L1R +
0L1R2L0L1L2L0R1R2L0R +
0L1R2L0L1L2L0R1R2R +
0L1R2L0L1L2R0L k
0L1R2L0L1L2R0R +
0L1R2L0L1R2L0L1L2L0L k
0L1R2L0L1R2L0L1L2L0R +
0L1R2L0L1R2L0L1L2R0L k
0L1R2L0L1R2L0L1L2R0R +
0L1R2L0L1R2L0L1R U
0L1R2L0L1R2L0R +
0L1R2L0L1R2R +
0L1R2L0R1L +
0L1R2L0R1R2L +
0L1R2L0R1R2R0L1L +
0L1R2L0R1R2R0L1R M
0L1R2L0R1R2R0R M
0L1R2R0L +
0L1R2R0R1L2L0L +
0L1R2R0R1L2L0R1L +
0L1R2R0R1L2L0R1R2L0L +
0L1R2R0R1L2L0R1R2L0R1L +
0L1R2R0R1L2L0R1R2L0R1R M
0L1R2R0R1L2L0R1R2R0L +
0L1R2R0R1L2L0R1R2R0R1L +
0L1R2R0R1L2L0R1R2R0R1R M
0L1R2R0R1L2R +
0L1R2R0R1R2L0L1L +
0L1R2R0R1R2L0L1R M
0L1R2R0R1R2L0R M
0L1R2R0R1R2R0L1L +
0L1R2R0R1R2R0L1R M
0L1R2R0R1R2R0R M
0R1L2L0L1L2L0L1L k
0R1L2L0L1L2L0L1R2L0L1L k
0R1L2L0L1L2L0L1R2L0L1R +
0R1L2L0L1L2L0L1R2L0R +
0R1L2L0L1L2L0L1R2R +
0R1L2L0L1L2L0R +
0R1L2L0L1L2R +
0R1L2L0L1R +
0R1L2L0R1L +
0R1L2L0R1R2L +
0R1L2L0R1R2R0L1L +
0R1L2L0R1R2R0L1R2L +
0R1L2L0R1R2R0L1R2R0L +
0R1L2L0R1R2R0L1R2R0R1L +
0R1L2L0R1R2R0L1R2R0R1R2L +
0R1L2L0R1R2R0L1R2R0R1R2R0L1L +
0R1L2L0R1R2R0L1R2R0R1R2R0L1R M
0R1L2L0R1R2R0L1R2R0R1R2R0R M
0R1L2L0R1R2R0R1L2L +
0R1L2L0R1R2R0R1L2R0L +
0R1L2L0R1R2R0R1L2R0R1L +
0R1L2L0R1R2R0R1L2R0R1R2L +
0R1L2L0R1R2R0R1L2R0R1R2R0L +
0R1L2L0R1R2R0R1L2R0R1R2R0R1L +
0R1L2L0R1R2R0R1L2R0R1R2R0R1R2L +
0R1L2L0R1R2R0R1L2R0R1R2R0R1R2R0L1L +
0R1L2L0R1R2R0R1L2R0R1R2R0R1R2R0L1R M
0R1L2L0R1R2R0R1L2R0R1R2R0R1R2R0R M
0R1L2L0R1R2R0R1R2L +
0R1L2L0R1R2R0R1R2R0L1L2L +
0R1L2L0R1R2R0R1R2R0L1L2R0L1L +
0R1L2L0R1R2R0R1R2R0L1L2R0L1R2L +
0R1L2L0R1R2R0R1R2R0L1L2R0L1R2R0L +
0R1L2L0R1R2R0R1R2R0L1L2R0L1R2R0R M
0R1L2L0R1R2R0R1R2R0L1L2R0R1L +
0R1L2L0R1R2R0R1R2R0L1L2R0R1R M
0R1L2L0R1R2R0R1R2R0L1R M
0R1L2L0R1R2R0R1R2R0R M
0R1L2R0L +
0R1L2R0R1L +
0R1L2R0R1R2L0L1L +
0R1L2R0R1R2L0L1R2L0L1L +
0R1L2R0R1R2L0L1R2L0L1R2L0L +
0R1L2R0R1R2L0L1R2L0L1R2L0R1L +
0R1L2R0R1R2L0L1R2L0L1R2L0R1R M
0R1L2R0R1R2L0L1R2L0L1R2R0L +
0R1L2R0R1R2L0L1R2L0L1R2R0R1L +
0R1L2R0R1R2L0L1R2L0L1R2R0R1R M
0R1L2R0R1R2L0L1R2L0R1L +
0R1L2R0R1R2L0L1R2L0R1R2L0L1L2L +
0R1L2R0R1R2L0L1R2L0R1R2L0L1L2R0L +
0R1L2R0R1R2L0L1R2L0R1R2L0L1L2R0R M
0R1L2R0R1R2L0L1R2L0R1R2L0L1R M
0R1L2R0R1R2L0L1R2L0R1R2L0R M
0R1L2R0R1R2L0L1R2L0R1R2R0L1L2L0L +
0R1L2R0R1R2L0L1R2L0R1R2R0L1L2L0R M
0R1L2R0R1R2L0L1R2L0R1R2R0L1L2R +
0R1L2R0R1R2L0L1R2L0R1R2R0L1R M
0R1L2R0R1R2L0L1R2L0R1R2R0R M
0R1L2R0R1R2L0L1R2R0L1L +
0R1L2R0R1R2L0L1R2R0L1R2L0L +
0R1L2R0R1R2L0L1R2R0L1R2L0R1L +
0R1L2R0R1R2L0L1R2R0L1R2L0R1R M
0R1L2R0R1R2L0L1R2R0L1R2R +
0R1L2R0R1R2L0L1R2R0R1L +
0R1L2R0R1R2L0L1R2R0R1R2L0L1L +
0R1L2R0R1R2L0L1R2R0R1R2L0L1R M
0R1L2R0R1R2L0L1R2R0R1R2L0R M
0R1L2R0R1R2L0L1R2R0R1R2R0L1L +
0R1L2R0R1R2L0L1R2R0R1R2R0L1R M
0R1L2R0R1R2L0L1R2R0R1R2R0R M
0R1L2R0R1R2L0R1L2L0L +
0R1L2R0R1R2L0R1L2L0R1L +
0R1L2R0R1R2L0R1L2L0R1R2L0L +
0R1L2R0R1R2L0R1L2L0R1R2L0R1L +
0R1L2R0R1R2L0R1L2L0R1R2L0R1R2L0L1L +
0R1L2R0R1R2L0R1L2L0R1R2L0R1R2L0L1R M
0R1L2R0R1R2L0R1L2L0R1R2L0R1R2L0R M
0R1L2R0R1R2L0R1L2L0R1R2L0R1R2R0L1L +
0R1L2R0R1R2L0R1L2L0R1R2L0R1R2R0L1R M
0R1L2R0R1R2L0R1L2L0R1R2L0R1R2R0R M
0R1L2R0R1R2L0R1L2L0R1R2R0L +
0R1L2R0R1R2L0R1L2L0R1R2R0R1L +
0R1L2R0R1R2L0R1L2L0R1R2R0R1R2L0L1L +
0R1L2R0R1R2L0R1L2L0R1R2R0R1R2L0L1R M
0R1L2R0R1R2L0R1L2L0R1R2R0R1R2L0R M
0R1L2R0R1R2L0R1L2L0R1R2R0R1R2R0L1L +
0R1L2R0R1R2L0R1L2L0R1R2R0R1R2R0L1R M
0R1L2R0R1R2L0R1L2L0R1R2R0R1R2R0R M
0R1L2R0R1R2L0R1L2R0L +
0R1L2R0R1R2L0R1L2R0R1L +
0R1L2R0R1R2L0R1L2R0R1R2L0L +
0R1L2R0R1R2L0R1L2R0R1R2L0R1L +
0R1L2R0R1R2L0R1L2R0R1R2L0R1R2L0L +
0R1L2R0R1R2L0R1L2R0R1R2L0R1R2L0R M
0R1L2R0R1R2L0R1L2R0R1R2L0R1R2R +
0R1L2R0R1R2L0R1L2R0R1R2R +
0R1L2R0R1R2L0R1R2L0L1L2L0L1L +
0R1L2R0R1R2L0R1R2L0L1L2L0L1R2L0L1L +
0R1L2R0R1R2L0R1R2L0L1L2L0L1R2L0L1R M
0R1L2R0R1R2L0R1R2L0L1L2L0L1R2L0R M
0R1L2R0R1R2L0R1R2L0L1L2L0L1R2R0L1L +
0R1L2R0R1R2L0R1R2L0L1L2L0L1R2R0L1R M
0R1L2R0R1R2L0R1R2L0L1L2L0L1R2R0R M
0R1L2R0R1R2L0R1R2L0L1L2L0R1L2L0L +
0R1L2R0R1R2L0R1R2L0L1L2L0R1L2L0R1L +
0R1L2R0R1R2L0R1R2L0L1L2L0R1L2L0R1R M
0R1L2R0R1R2L0R1R2L0L1L2L0R1L2R0L +
0R1L2R0R1R2L0R1R2L0L1L2L0R1L2R0R1L +
0R1L2R0R1R2L0R1R2L0L1L2L0R1L2R0R1R M
0R1L2R0R1R2L0R1R2L0L1L2L0R1R M
0R1L2R0R1R2L0R1R2L0L1L2R0L1L +
0R1L2R0R1R2L0R1R2L0L1L2R0L1R2L0L1L +
0R1L2R0R1R2L0R1R2L0L1L2R0L1R2L0L1R M
0R1L2R0R1R2L0R1R2L0L1L2R0L1R2L0R M
0R1L2R0R1R2L0R1R2L0L1L2R0L1R2R0L1L +
0R1L2R0R1R2L0R1R2L0L1L2R0L1R2R0L1R M
0R1L2R0R1R2L0R1R2L0L1L2R0L1R2R0R M
0R1L2R0R1R2L0R1R2L0L1L2R0R1L2L0L +
0R1L2R0R1R2L0R1R2L0L1L2R0R1L2L0R1L +
0R1L2R0R1R2L0R1R2L0L1L2R0R1L2L0R1R M
0R1L2R0R1R2L0R1R2L0L1L2R0R1L2R0L +
0R1L2R0R1R2L0R1R2L0L1L2R0R1L2R0R1L +
0R1L2R0R1R2L0R1R2L0L1L2R0R1L2R0R1R M
0R1L2R0R1R2L0R1R2L0L1L2R0R1R M
0R1L2R0R1R2L0R1R2L0L1R M
0R1L2R0R1R2L0R1R2L0R M
0R1L2R0R1R2L0R1R2R0L1L2L0L1L +
0R1L2R0R1R2L0R1R2R0L1L2L0L1R2L0L +
0R1L2R0R1R2L0R1R2R0L1L2L0L1R2L0R M
0R1L2R0R1R2L0R1R2R0L1L2L0L1R2R +
0R1L2R0R1R2L0R1R2R0L1L2L0R1L +
0R1L2R0R1R2L0R1R2R0L1L2L0R1R M
0R1L2R0R1R2L0R1R2R0L1L2R0L +
0R1L2R0R1R2L0R1R2R0L1L2R0R1L +
0R1L2R0R1R2L0R1R2R0L1L2R0R1R M
0R1L2R0R1R2L0R1R2R0L1R M
0R1L2R0R1R2L0R1R2R0R M
0R1L2R0R1R2R0L +
0R1L2R0R1R2R0R1L +
0R1L2R0R1R2R0R1R2L0L1L +
0R1L2R0R1R2R0R1R2L0L1R M
0R1L2R0R1R2R0R1R2L0R M
0R1L2R0R1R2R0R1R2R +
0R1R2L0L1L2L +
0R1R2L0L1L2R0L1L +
0R1R2L0L1L2R0L1R M
0R1R2L0L1L2R0R1L2L +
0R1R2L0L1L2R0R1L2R0L1L +
0R1R2L0L1L2R0R1L2R0L1R M
0R1R2L0L1L2R0R1L2R0R1L2L +
0R1R2L0L1L2R0R1L2R0R1L2R0L +
0R1R2L0L1L2R0R1L2R0R1L2R0R M
0R1R2L0L1L2R0R1L2R0R1R M
0R1R2L0L1L2R0R1R M
0R1R2L0L1R M
0R1R2L0R M
0R1R2R0L1L2L0L1L +
0R1R2R0L1L2L0L1R M
0R1R2R0L1L2L0R1L2L0L1L +
0R1R2R0L1L2L0R1L2L0L1R M
0R1R2R0L1L2L0R1L2L0R1L2L0L1L +
0R1R2R0L1L2L0R1L2L0R1L2L0L1R M
0R1R2R0L1L2L0R1L2L0R1L2L0R M
0R1R2R0L1L2L0R1L2L0R1L2R0L1L +
0R1R2R0L1L2L0R1L2L0R1L2R0L1R M
0R1R2R0L1L2L0R1L2L0R1L2R0R M
0R1R2R0L1L2L0R1L2L0R1R M
0R1R2R0L1L2L0R1L2R0L1L +
0R1R2R0L1L2L0R1L2R0L1R M
0R1R2R0L1L2L0R1L2R0R1L2L0L1L +
0R1R2R0L1L2L0R1L2R0R1L2L0L1R M
0R1R2R0L1L2L0R1L2R0R1L2L0R M
0R1R2R0L1L2L0R1L2R0R1L2R0L +
0R1R2R0L1L2L0R1L2R0R1L2R0R M
0R1R2R0L1L2L0R1L2R0R1R M
0R1R2R0L1L2L0R1R M
0R1R2R0L1L2R0L1L +
0R1R2R0L1L2R0L1R M
0R1R2R0L1L2R0R1L2L0L1L +
0R1R2R0L1L2R0R1L2L0L1R M
0R1R2R0L1L2R0R1L2L0R1L +
0R1R2R0L1L2R0R1L2L0R1R M
0R1R2R0L1L2R0R1L2R +
0R1R2R0L1L2R0R1R M
0R1R2R0L1R M
0R1R2R0R M"""


def partition():
    trie = {}
    parsed = []
    for row in LEAVES.splitlines():
        path, label = row.split()
        tokens = [(int(path[j]), path[j+1]) for j in range(0, len(path), 2)]
        node = trie
        for token in tokens:
            assert 'leaf' not in node
            node = node.setdefault(token, {})
        assert not node
        node['leaf'] = label
        parsed.append((tokens, label))

    def full(node):
        if 'leaf' in node:
            assert len(node) == 1
            return
        assert len(node) == 2
        axes = {token[0] for token in node}
        assert len(axes) == 1
        axis = next(iter(axes))
        assert set(node) == {(axis, 'L'), (axis, 'R')}
        for child in node.values():
            full(child)
    full(trie)
    assert len(parsed) == 238
    assert max(len(path) for path, _ in parsed) == 17
    return parsed


def main():
    qs = exact_coefficients(9)
    ell = [None] + [add(scale(A, Q(1, k+1)), scale(C, Q(1, k)))
                    for k in range(1, 5)]
    assert qs[1] == ell[1]
    assert qs[2] == add(ell[2], scale(power(ell[1], 2), Q(1, 2)))
    assert qs[3] == add(ell[3], mul(ell[1], ell[2]),
                        scale(power(ell[1], 3), Q(1, 6)))
    assert qs[4] == add(ell[4], mul(ell[1], ell[3]),
                        scale(power(ell[2], 2), Q(1, 2)),
                        scale(mul(power(ell[1], 2), ell[2]), Q(1, 2)),
                        scale(power(ell[1], 4), Q(1, 24)))
    ds = [add(*(scale(qs[j+1], Q((-1)**j * comb(n, j)))
                 for j in range(n+1))) for n in range(9)]
    assert ds[0] == qs[1] and ds[1] == add(qs[1], scale(qs[2], -1))
    print('FIDELITY: source q1,...,q4 identities and q0,...,q9 generating recurrence: OK')
    print('FIDELITY: source binomial d0=q1, d1=q1-q2 and d0,...,d8 reconstruction: OK')

    # P6 = a^2 sum K^(i)(0)t^i/i!, without stored polynomial coefficients.
    kminus2 = add(A, C, scale(ONE, -2))
    kminus1 = add(A, C, scale(ONE, -1))
    middle = mul(A, add(scale(A, 2), scale(C, 2), scale(ONE, -3)))
    P6 = {}
    for i in range(7):
        pi = add(mul(power(A, 2), ds[i+2]), scale(mul(middle, ds[i+1]), -1),
                 mul(mul(kminus2, kminus1), ds[i]))
        P6 = add(P6, {(x, y, z+i): v / factorial(i) for (x,y,z),v in pi.items()})
    assert tuple(max(e[i] for e in P6) for i in range(3)) == (11, 9, 6)

    # Rigorous logarithmic relaxation of delta = 3 log(3)-2.
    delta0 = Q(1295836866004, 10**12)
    log_lo = 2 * sum((Q(1, (2*j+1) * 2**(2*j+1)) for j in range(24)), Q(0))
    log_tail = Q(2, 49 * 2**49) / (1-Q(1, 4))
    assert delta0-Q(1,10**12) < 3*log_lo-2
    assert 3*(log_lo+log_tail)-2 < delta0+Q(1,10**12)
    assert 1 < 3*log_lo-2 < 3*(log_lo+log_tail)-2 < Q(13,10)
    derivative_bound = 2*Q(7,3)*(5*Q(7,3)+26) + 2*Q(7,3)*(Q(13,10)*Q(7,3)+4)*(Q(7,3)+4)
    assert derivative_bound == Q(51793,135) < 500
    assert 500*Q(1,10**12) < Q(1,10**9)
    ac4 = add(A, scale(C, 4))
    v = add(scale(A, delta0), scale(C, 4))
    Mbar = add(scale(mul(v, add(scale(A, 5), scale(C, 26))), 2),
               scale(power(ac4, 2), -8), scale(mul(power(v, 2), ac4), -1),
               scale(ONE, Q(1,10**9)))
    U = add(scale(add(ONE, scale(C,-1)), 4),
            scale(power(add(A, scale(ONE,-2), scale(C,2)),2),-1))
    print('LOG_RELAXATION: |delta-delta0|<1/10^12 and derivative bound 51793/135<500: OK')

    root = ((Q(1),Q(7,3)), (Q(0),Q(1)), (Q(0),Q(16)))
    rho = Q(6583,11943936)
    minimum = None
    counts = Counter()
    for tokens, label in partition():
        box = list(root)
        for axis, side in tokens:
            lo, hi = box[axis]
            mid = (lo+hi)/2
            box[axis] = (lo, mid) if side=='L' else (mid, hi)
        poly = {'+':P6, 'M':Mbar, 'U':U, 'k':kminus2}[label]
        values = bernstein_values(poly, tuple(box))
        if label=='+':
            lower = min(values)
            assert lower >= rho, (tokens, lower)
            minimum = lower if minimum is None else min(minimum,lower)
        else:
            assert max(values) < 0, (tokens,label,max(values))
        counts[label] += 1
    assert counts == {'+':134,'M':87,'U':1,'k':16}
    assert minimum == rho
    print('PARTITION: 238 closed leaves; maximum depth 17; 134 positive, 87 M, 1 U, 16 k: OK')
    print('UNIFORM_FINITE_INTERVAL: minimum Bernstein coefficient =', minimum)

    tail = {}
    for (x,y,z),v in P6.items():
        if z:
            e=(x,y,0)
            tail[e]=tail.get(e,Q(0))+v*z*16**(z-1)
    tail={e:v for e,v in tail.items() if v}
    eta = min(bernstein_values(tail, root))
    assert eta == Q(45599851,571536000) > 0
    print('TAIL: minimum Bernstein coefficient of dP6/dt at t=16 =', eta)
    print('SCOPE: exact polynomial premises verified; source integral and analytic bridges require Lean proof')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
