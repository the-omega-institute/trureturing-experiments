#!/usr/bin/env python3
"""Exact graph and derivative checks for Dong et al., Conjecture 2.

Python 3 with SymPy. No floating-point arithmetic or numerical roots are used.
The connectivity transfer computes Whitney's spanning-subgraph polynomial from
an actual edge list, independently of the endpoint-coloring closed formula.
"""
from fractions import Fraction
from itertools import combinations, product
from math import factorial
import sympy as sp

q, t = sp.symbols('q t')


def canonical(labels):
    renaming = {}
    return tuple(renaming.setdefault(x, len(renaming)) for x in labels)


def graph_polynomial(n, edges):
    """Sum (-1)^|A| q^components(V,A), with frontier partitions.

    Each edge is absent or present. A present edge merges its blocks and
    reverses its weight's sign. Forgetting a completed component supplies q.
    """
    assert len(set(edges)) == len(edges)
    assert all(0 <= u < v < n for u, v in edges)
    last = {v: i for i, e in enumerate(edges) for v in e}
    active, states = [], {(): sp.Integer(1)}
    for i, (u, v) in enumerate(edges):
        extended = active + [x for x in (u, v) if x not in active]
        keep = [j for j, x in enumerate(extended) if last[x] > i]
        new_states = {}
        for partition, weight in states.items():
            labels = list(partition)
            for _ in range(len(extended) - len(active)):
                labels.append(max(labels, default=-1) + 1)
            for present in (False, True):
                lab = labels.copy()
                if present:
                    old, new = lab[extended.index(v)], lab[extended.index(u)]
                    lab = [new if a == old else a for a in lab]
                closed = len(set(lab) - {lab[j] for j in keep})
                key = canonical([lab[j] for j in keep])
                contribution = (-1 if present else 1) * weight * q**closed
                new_states[key] = new_states.get(key, 0) + contribution
        states = {s: sp.expand(w) for s, w in new_states.items()}
        active = [extended[j] for j in keep]
    assert not active and set(states) == {()}
    return sp.Poly(states[()] * q**(n - len(last)), q, domain=sp.ZZ).as_expr()


def theta(a, copies=1):
    """Three internally disjoint a-edge paths; copies share only vertex 0."""
    assert a >= 2 and copies >= 1
    edges, next_vertex = [], 1
    for _ in range(copies):
        endpoint = next_vertex
        next_vertex += 1
        for _ in range(3):
            path = [0] + list(range(next_vertex, next_vertex + a - 1)) + [endpoint]
            next_vertex += a - 1
            edges.extend(tuple(sorted(e)) for e in zip(path, path[1:]))
    return next_vertex, edges


def connected(n, edges):
    seen = {0}
    while True:
        extended = seen | {v for u, v in edges if u in seen} | {u for u, v in edges if v in seen}
        if extended == seen:
            return len(seen) == n
        seen = extended


def brute_colorings(n, edges, colors):
    return sum(all(f[u] != f[v] for u, v in edges)
               for f in product(range(colors), repeat=n))


def log_jet(poly, point, degree):
    """Derivatives of log(poly), by the exact identity P'=P (log P)'."""
    coeff = [Fraction(sp.diff(poly, t, j).subs(t, point)) / factorial(j)
             for j in range(degree + 1)]
    assert coeff[0] > 0
    log_coeff = [Fraction(0)] * (degree + 1)
    for j in range(1, degree + 1):
        residue = sum(coeff[i] * (j-i) * log_coeff[j-i] for i in range(1, j))
        log_coeff[j] = (j * coeff[j] - residue) / (j * coeff[0])
    return [factorial(j) * log_coeff[j] for j in range(degree + 1)]


def main():
    # Independent coloring-definition cross-check for every graph on 1..4 vertices.
    graph_count = 0
    for n in range(1, 5):
        possible = list(combinations(range(n), 2))
        for mask in range(1 << len(possible)):
            edges = [e for i, e in enumerate(possible) if mask >> i & 1]
            P = graph_polynomial(n, edges)
            for colors in range(1, 5):
                assert P.subs(q, colors) == brute_colorings(n, edges, colors)
            graph_count += 1
    print(f'COLORING_FIDELITY: {graph_count} graphs, colors=1..4')

    # A source-stated result: Section 1, equation (7), epsilon(K_n)=n-H_n.
    # The preceding proposition gives epsilon(G)=n+P\'(G,-1)/P(G,-1).
    for n in range(1, 7):
        P = graph_polynomial(n, list(combinations(range(n), 2)))
        assert sp.expand(P - sp.prod(q-i for i in range(n))) == 0
        epsilon = n + Fraction(sp.diff(P, q).subs(q, -1), P.subs(q, -1))
        assert epsilon == n - sum((Fraction(1, i) for i in range(1, n+1)), Fraction(0))
    print('SOURCE_FIDELITY: epsilon(K_n)=n-H_n, n=1..6 (Section 1, equation (7))')

    polynomials = {}
    for a in range(2, 10):
        n, edges = theta(a)
        assert n == 3*a-1 and len(edges) == 3*a and connected(n, edges)
        P = graph_polynomial(n, edges)
        S = sp.cancel(((q-1)**a + (q-1)*(-1)**a)/q)
        D = sp.cancel(((q-1)**a - (-1)**a)/q)
        assert sp.cancel(P - (q*S**3 + q*(q-1)*D**3)) == 0
        Q = sp.expand((-1)**n * P.subs(q, -t))
        closed = (1+t)*((1+t)**(3*a-1)-3*(1+t)**a+2+t)/t
        assert sp.cancel(Q-closed) == 0
        assert sp.degree(Q, t) == n
        assert all(c >= 0 for c in sp.Poly(Q, t).all_coeffs())
        polynomials[a] = P, Q
    print('THETA_FIDELITY: actual graph polynomial = endpoint formula = Q_a, a=2..9')

    P, Q = polynomials[7]
    jets = tuple(int(sp.diff(Q, t, j).subs(t, 1)) for j in range(6))
    expected_jets = (2096390, 19920639, 180348932, 1550830068, 12623781168, 96911980560)
    assert jets == expected_jets
    p, a, b, c, d, e = jets
    fifth_formula = Fraction(-e*p**4 + (5*d*a+10*c*b)*p**3
                            - (20*c*a*a+30*b*b*a)*p*p + 60*b*a**3*p - 24*a**5, p**5)
    fifth_recurrence = -log_jet(Q, 1, 5)[5]  # x=-t contributes (-1)^5.
    # Third route: differentiate Q'/Q four times, using exact symbolic algebra.
    fifth_symbolic = -sp.diff(sp.diff(Q, t)/Q, t, 4).subs(t, 1)
    expected = Fraction(558835859111752710016430691303, 5061397078105038614084902487500)
    assert fifth_formula == fifth_recurrence == fifth_symbolic == expected > 0
    print('Q_7_JETS_AT_1:', jets)
    print('L_THETA7_FIFTH_AT_MINUS1:', fifth_formula)
    print('DERIVATIVE_METHODS: closed formula, rational Taylor recurrence, symbolic differentiation agree')

    # Reconstruct the actual glued graphs independently for several copy counts.
    for m in range(1, 7):
        n, edges = theta(7, m)
        assert n == 19*m+1 and len(edges) == 21*m and connected(n, edges)
        Pm = graph_polynomial(n, edges)
        assert sp.Poly(Pm, q) * sp.Poly(q**(m-1), q) == sp.Poly(P, q)**m
        Qm = sp.expand((-1)**n * Pm.subs(q, -t))
        value = -log_jet(Qm, 1, 5)[5]
        assert value == m*expected + 24*(m-1) > 0
        print(f'GLUED_FAMILY: m={m}, vertices={n}, edges={len(edges)}, fifth={value}')
    # Exact symbolic positivity certificate for every integer m>=1:
    # m*c+24(m-1)=c+(m-1)(c+24), where both coefficients are positive.
    assert expected > 0 and expected + 24 > 0
    print('ALL_M_CERTIFICATE: c+(m-1)*(c+24)>0 for every integer m>=1')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
