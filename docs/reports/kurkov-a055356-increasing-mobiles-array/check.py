#!/usr/bin/env python3
"""Exact finite checks for Kurkov's final FORMULA conjecture in OEIS A055356.

All coefficients and combinatorial weights are Python integers. This program
checks finite cases; it does not certify the universal conjecture.
"""
from collections import defaultdict
from functools import cache
from itertools import product
from math import comb, factorial


def trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return tuple(p)


def add(p, q):
    return trim([(p[i] if i < len(p) else 0) +
                 (q[i] if i < len(q) else 0) for i in range(max(len(p), len(q)))])


def scale(c, p):
    return trim([c * v for v in p])


def mul(p, q):
    out = [0] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i + j] += a * b
    return trim(out)


def xmul(p):
    return trim([0, *p])


def xderiv(p):
    return trim([i * a for i, a in enumerate(p)])


@cache
def T(N, leaves):
    if N < 1 or leaves < 1 or leaves > N:
        return 0
    if N == 1:
        return int(leaves == 1)
    return (N - 2) * T(N - 1, leaves - 1) + leaves * T(N - 1, leaves)


def P(N):
    return trim([T(N, leaves) for leaves in range(1, N + 1)])


@cache
def R(n, k):
    if n == 0:
        return (1,)
    convolution = (0,)
    for i in range(n):
        for j in range(k + 1):
            term = scale(comb(n - 1, i), mul(R(n - i - 1, j), R(i, k - j)))
            convolution = add(convolution, term)
    return add(R(n - 1, k + 1), xmul(convolution))


def derive(poly):
    """D on Z[x,a,b]: D(x)=0, D(a)=ab+xa^2, D(b)=xab.

    Keys are exponents (x,a,b). No evaluated-row recurrence is used here.
    """
    out = defaultdict(int)
    for (r, s, t), c in poly.items():
        if s:
            out[r, s, t + 1] += c * s
            out[r + 1, s + 1, t] += c * s
        if t:
            out[r + 1, s + 1, t] += c * t
    return {m: c for m, c in out.items() if c}


def evaluate(poly):
    out = [0] * (1 + max(m[0] for m in poly))
    for (r, _s, _t), c in poly.items():
        out[r] += c
    return trim(out)


def linear_next(n, k, p):
    return add(add(p, scale(n + k + 1, xmul(p))), xderiv(p))


def mobile_counts(N):
    """Enumerate increasing parent trees, weighting all cyclic child orders.

    For a vertex of outdegree d>0 there are (d-1)! cyclic child orders;
    the empty child set has one order. Labels are 0,...,N-1 with root 0.
    """
    counts = [0] * N
    for parents in product(*(range(v) for v in range(1, N))):
        degree = [0] * N
        for parent in parents:
            degree[parent] += 1
        leaves = sum(d == 0 for d in degree[1:])
        weight = 1
        for d in degree:
            if d:
                weight *= factorial(d - 1)
        counts[leaves - 1] += weight
    return tuple(counts)


def ma_grammar_third_derivative():
    """Ma's G={u->u^2 v,v->uv}; keys are (u,v) exponents."""
    p = {(0, 1): 1}
    for _ in range(3):
        out = defaultdict(int)
        for (s, t), c in p.items():
            if s:
                out[s + 1, t + 1] += c * s
            if t:
                out[s + 1, t] += c * t
        p = dict(out)
    return p


def main():
    # Literal OEIS example rows, including the terminal zeros.
    published = [(1,), (1, 0), (1, 1, 0), (1, 4, 2, 0),
                 (1, 11, 18, 6, 0), (1, 26, 98, 96, 24, 0),
                 (1, 57, 424, 874, 600, 120, 0)]
    for N, row in enumerate(published, 1):
        assert tuple(T(N, k) for k in range(1, N + 1)) == row
    print('Fidelity: all 7 displayed OEIS rows, with terminal zeros, reproduced.')
    for N in range(2, 9):
        actual = mobile_counts(N)
        assert actual == tuple(T(N, k) for k in range(1, N + 1))
    print('Fidelity: independent parent-tree/cyclic-order counts agree for 2 <= N <= 8.')
    assert ma_grammar_third_derivative() == {(3, 1): 1, (3, 2): 4, (3, 3): 2}
    print('Indexing: Ma grammar D^3(v) = u^3*(v + 4*v^2 + 2*v^3), matching OEIS row N=4.')
    for N in range(1, 13):
        p = P(N)
        whitmore = add(add(p, scale(N - 1, xmul(p))), xderiv(p))
        assert P(N + 1) == whitmore
    print('Fidelity: Whitmore differential recurrence agrees for 1 <= N <= 12.')
    checked = 0
    for k in range(9):
        poly = {(0, 1, k): 1}
        linear = (1,)
        for n in range(11):
            for r, s, t in poly:
                assert 0 <= r <= n and s == r + 1 and t == n + k - r
            q = evaluate(poly)
            assert q == linear == R(n, k), (n, k, q, linear, R(n, k))
            if k == 0:
                assert q == P(n + 2), (n, q, P(n + 2))
            poly = derive(poly)
            linear = linear_next(n, k, linear)
            checked += 1
    print(f'Array: original nonlinear recurrence, direct derivation and linear differential recurrence agree on {checked} cases (0 <= n <= 10, 0 <= k <= 8).')
    print('Support: every direct derivation monomial has exponents (r,r+1,n+k-r), 0 <= r <= n.')
    print('Claim: P_(n+2) = R(n,0) for every checked n from 0 through 10.')
    for n in range(11):
        print(f'R({n},0) coefficients: {list(R(n, 0))}')
    print('Scope: finite exact checks only; no universal proof or Lean kernel check is claimed.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
