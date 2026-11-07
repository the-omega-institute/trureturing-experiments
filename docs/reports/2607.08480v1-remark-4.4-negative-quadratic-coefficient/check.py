#!/usr/bin/env python3
"""Exact checks for arXiv:2607.08480v1, Remark 4.4; no Lean proof.

Source fidelity: the source's star polynomial and Theorem 4.1 on all
labelled simple graphs of orders 1 through 5. The target witness is counted
from explicit adjacency, independently of the symbolic product identity.
"""
from functools import cache
from itertools import combinations, product
from math import comb

import sympy as sp


def add(a, b):
    out = [0] * max(len(a), len(b))
    for i, v in enumerate(a):
        out[i] += v
    for i, v in enumerate(b):
        out[i] += v
    return tuple(out)


def multiply(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return tuple(out)


def eval_at_minus_one(p):
    return sum((-1) ** i * v for i, v in enumerate(p))


def divide_one_plus_x(p):
    assert len(p) > 1 and eval_at_minus_one(p) == 0
    q = [p[0]]
    for c in p[1:-1]:
        q.append(c - q[-1])
    assert q[-1] == p[-1] and multiply((1, 1), q) == p
    return tuple(q)


def exact_multiplicity(p):
    assert p and p[0] == 1
    m = 0
    while eval_at_minus_one(p) == 0:
        p = divide_one_plus_x(p)
        m += 1
    return m


def adjacency(n, edges):
    a = [0] * n
    for u, v in edges:
        assert 0 <= u < v < n
        a[u] |= 1 << v
        a[v] |= 1 << u
    return tuple(a)


def independent_sequence(a):
    """Exhaustive deletion/contraction count on actual vertex subsets."""
    @cache
    def count(mask):
        if not mask:
            return (1,)
        bit = mask & -mask
        v = bit.bit_length() - 1
        rest = mask ^ bit
        absent = count(rest)
        present = (0,) + count(rest & ~a[v])
        return add(absent, present)
    return count((1 << len(a)) - 1)


def connected(a):
    reached = 1
    frontier = 1
    while frontier:
        bit = frontier & -frontier
        frontier ^= bit
        v = bit.bit_length() - 1
        new = a[v] & ~reached
        reached |= new
        frontier |= new
    return reached == (1 << len(a)) - 1


def clique_join_graph(sizes, k):
    blocks = [i for i, size in enumerate(sizes) for _ in range(size)]
    blocks += [len(sizes)] * k
    n = len(blocks)
    universal = len(sizes)
    edges = [(i, j) for i, j in combinations(range(n), 2)
             if blocks[i] == blocks[j] or universal in (blocks[i], blocks[j])]
    return adjacency(n, edges), len(edges)


def main():
    # Source fidelity, §2: P_{S_n}=(1+x)^(n-1)+x for n>=2.
    for n in range(2, 11):
        a = adjacency(n, [(i, n - 1) for i in range(n - 1)])
        expected = list(comb(n - 1, i) for i in range(n))
        expected[1] += 1
        assert independent_sequence(a) == tuple(expected)
    print('FIDELITY: star polynomial reproduced for 2 <= n <= 10')

    # Source fidelity, Theorem 4.1. This checks the established alpha-2 bound,
    # not the open alpha-3 positivity question, on every graph in this range.
    graphs = eligible = 0
    for n in range(1, 6):
        possible = tuple(combinations(range(n), 2))
        for mask in range(1 << len(possible)):
            a = adjacency(n, [e for i, e in enumerate(possible) if mask & (1 << i)])
            p = independent_sequence(a)
            alpha = len(p) - 1
            m = exact_multiplicity(p)
            graphs += 1
            if m >= alpha - 2:
                eligible += 1
                assert all(p[i] ** 2 >= p[i - 1] * p[i + 1]
                           for i in range(1, alpha))
    print(f'FIDELITY: Theorem 4.1 verified on {graphs} graphs; {eligible} meet its hypothesis')

    # Construct every edge of the t=5 graph explicitly.
    sizes, k = (3, 4, 5, 5), 96
    a, edges = clique_join_graph(sizes, k)
    raw = independent_sequence(a)
    p_product = (1,)
    for size in sizes:
        p_product = multiply(p_product, (1, size))
    p_product = add(p_product, (0, k))
    assert raw == p_product == (1, 113, 107, 295, 300)

    # A second brute-force enumeration lists every allowed clique choice.
    offsets, offset = [], 0
    for size in sizes:
        offsets.append(tuple(range(offset, offset + size)))
        offset += size
    counts = [0] * 5
    for choice in product(*[(None,) + block for block in offsets]):
        vertices = tuple(v for v in choice if v is not None)
        assert all(not (a[u] & (1 << v)) for u, v in combinations(vertices, 2))
        counts[len(vertices)] += 1
    counts[1] += k  # Only singletons in the joining clique can occur.
    assert tuple(counts) == raw
    assert len(a) == 113 and edges == 6221 and len(raw) - 1 == 4
    assert connected(a) and all(a)
    q = divide_one_plus_x(raw)
    assert q == (1, 112, -5, 300)
    assert eval_at_minus_one(q) == -416 and exact_multiplicity(raw) == 1
    assert comb(113 - 4 + 3, 2) - edges == q[2] == -5
    assert raw[1] > raw[2] < raw[3]
    print(f'WITNESS: n={len(a)}, edges={edges}, alpha=4, connected=True, no_isolates=True')
    print(f'WITNESS: independence_sequence={raw}; all {sum(raw)} independent sets counted')
    print(f'WITNESS: quotient={q}; Q(-1)={eval_at_minus_one(q)}; exact_multiplicity=1; b=-5')
    print('TRIAGE ONLY: strict coefficient valley 113 > 107 < 295')

    # Classical partial prior art has the same four blocks but joining size 25.
    old = add(multiply(multiply((1, 3), (1, 4)),
                       multiply((1, 5), (1, 5))), (0, 25))
    assert old == (1, 42, 107, 295, 300) and eval_at_minus_one(old) == 71
    print('PRIOR ART: Levit-Mandrescu joining size 25 gives P(-1)=71, not the required root')

    # Exact symbolic identities establish the family's algebraic certificates.
    x, t, u = sp.symbols('x t u', integer=True)
    h = (1 + 3*x) * (1 + 4*x) * (1 + t*x)**2
    k_t = 6 * (t - 1)**2
    p_t = sp.expand(h + k_t*x)
    n_t = 6*t**2 - 10*t + 13
    q_t = 1 + (n_t - 1)*x + t*(24 - 5*t)*x**2 + 12*t**2*x**3
    assert sp.expand(h.subs(x, -1) - k_t) == 0
    assert sp.expand(p_t - (1 + x)*q_t) == 0
    assert sp.expand(q_t.subs(x, -1) + (t - 1)*(23*t - 11)) == 0
    assert [sp.expand(p_t).coeff(x, i) for i in range(5)] == [
        1, n_t, t**2 + 14*t + 12, 7*t**2 + 24*t, 12*t**2]
    # Edge count from the five cliques plus every cross edge of the join.
    e_t = 3 + 6 + t*(t - 1) + k_t*(k_t - 1)/sp.Integer(2) + k_t*(7 + 2*t)
    source_b = (n_t - 1)*(n_t - 2)/sp.Integer(2) - e_t
    source_c = e_t + ((n_t-1)*(n_t-2)*(n_t-3) - n_t*(n_t-1)*(n_t-2))/sp.Integer(6) + 7*t**2 + 24*t
    assert sp.expand(source_b - t*(24 - 5*t)) == 0
    assert sp.expand(source_c - 12*t**2) == 0
    gap1 = sp.expand(p_t.coeff(x, 1) - p_t.coeff(x, 2))
    gap3 = sp.expand(p_t.coeff(x, 3) - p_t.coeff(x, 2))
    certificates = {
        '-b': sp.expand((-t*(24 - 5*t)).subs(t, u + 5)),
        '-Q(-1)': sp.expand(((t - 1)*(23*t - 11)).subs(t, u + 5)),
        's1-s2': sp.expand(gap1.subs(t, u + 5)),
        's3-s2': sp.expand(gap3.subs(t, u + 5)),
    }
    for name, polynomial in certificates.items():
        coeffs = sp.Poly(polynomial, u).all_coeffs()
        assert all(c > 0 for c in coeffs)
        print(f'FAMILY: {name} at t=u+5 is {polynomial}, positive for every integer u>=0')
    print('FAMILY: H(-1)=6(t-1)^2; P=(1+x)Q; Q(-1)=-(t-1)(23t-11); b=t(24-5t)')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
