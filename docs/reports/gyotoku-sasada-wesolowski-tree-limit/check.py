#!/usr/bin/env python3
"""Exact finite sanity checks for arXiv:2609.15199v1; no numerical limit evidence.

Uses SymPy exact polynomial/rational algebra. The finite checks do not prove the
universal analytic limit, integrability, or novelty of the proposed settlement.
"""
from itertools import product
from math import prod
import sympy as s


def trees(n):
    for word in product(range(n), repeat=n - 2):
        degree = [1] * n
        for v in word:
            degree[v] += 1
        edges = []
        for v in word:
            leaf = next(i for i in range(n) if degree[i] == 1)
            edges.append(tuple(sorted((leaf, v))))
            degree[leaf] -= 1
            degree[v] -= 1
        a, b = [i for i in range(n) if degree[i] == 1]
        yield tuple(sorted(edges + [(a, b)]))


def root_tree(n, edges, root):
    adj = [[] for _ in range(n)]
    for a, b in edges:
        adj[a].append(b)
        adj[b].append(a)
    parent, order = {root: None}, [root]
    for v in order:
        for w in adj[v]:
            if w not in parent:
                parent[w] = v
                order.append(w)
    child = {v: [w for w in adj[v] if parent.get(w) == v] for v in range(n)}
    return child, order


def delta(n, edges, q, k, vertices=None):
    vertices = list(range(n)) if vertices is None else list(vertices)
    terms = []
    for mask in range(1 << len(vertices)):
        u = {vertices[i] for i in range(len(vertices)) if mask & (1 << i)}
        terms.append((-1) ** len(u) * prod(k[v] for v in u)
                     * prod(q[e] for e in edges if set(e) <= u))
    return s.expand(sum(terms))


def zero(expr):
    assert s.cancel(expr) == 0, expr


def check():
    eps = s.Symbol('eps', positive=True)
    cases = roots = minors = 0
    for n in range(2, 5):
        x = s.symbols(f'x0:{n}')
        y = s.symbols(f'y0:{n}')
        z = s.symbols(f'z0:{n}', positive=True)
        b = s.symbols(f'b0:{n}', positive=True)
        for edges in trees(n):
            cases += 1
            c = dict(zip(edges, s.symbols(f'c0:{n-1}', nonzero=True, real=True)))
            q = dict(zip(edges, s.symbols(f'q0:{n-1}')))
            qe = {e: 1 - eps**2 * c[e]**2 for e in edges}
            # All induced principal polynomials, including the empty one.
            for mask in range(1 << n):
                u = [v for v in range(n) if mask & (1 << v)]
                matrix = s.Matrix([[y[v] if v == w else c.get(tuple(sorted((v, w))), 0)
                                    for w in u] for v in u])
                polynomial = s.Poly(delta(n, edges, qe, [1-eps*t for t in y], u), eps)
                for j in range(len(u)):
                    zero(polynomial.nth(j))
                zero(polynomial.nth(len(u)) - matrix.det())
                minors += 1
            for root in range(n):
                roots += 1
                child, order = root_tree(n, edges, root)
                # Fidelity: Lemma 2.5 inverse, Lemma 2.3 recursion, Lemma 2.4 factorization.
                phi = [x[v] * prod((1-x[w])/(1-q[tuple(sorted((v, w)))]*x[w])
                                  for w in child[v]) for v in range(n)]
                psi = {}
                for v in reversed(order):
                    psi[v] = s.cancel(phi[v] * prod(
                        (1-q[tuple(sorted((v, w)))]*psi[w])/(1-psi[w]) for w in child[v]))
                    zero(psi[v] - x[v])
                zero(delta(n, edges, q, phi) - prod(1-t for t in x))
                # Expected settlement: exact first-order inverse transport and limiting Schur pivots.
                h = [z[v] + sum(c[tuple(sorted((v, w)))]**2/z[w] for w in child[v])
                     for v in range(n)]
                matrix = s.Matrix([[h[v] if v == w else c.get(tuple(sorted((v, w))), 0)
                                    for w in range(n)] for v in range(n)])
                zero(matrix.det() - prod(z))
                zero(s.Matrix(h).jacobian(z).det() - 1)
                zero(sum(b[v]*h[v] for v in range(n)) - sum(b[v]*z[v] for v in range(n))
                     - sum(b[v]*c[tuple(sorted((v,w)))]**2/z[w]
                           for v in range(n) for w in child[v]))
                for v in range(n):
                    rv = prod(z[w]/(z[w]+eps*c[tuple(sorted((v,w)))]**2*(1-eps*z[w]))
                              for w in child[v])
                    he = s.cancel((1-(1-eps*z[v])*rv)/eps)
                    zero(he.subs(eps, 0) - h[v])
                    # Exact diagonal Jacobian: this extra factor is essential for density L1.
                    zero(s.diff(he, z[v]) - rv)
    print(f'Exact source fidelity: {roots} rooted labeled trees, all labeled trees of sizes 2..4')
    print(f'Exact scaling/minor checks: {minors} induced minors on {cases} labeled trees')
    print(f'Exact MY determinant, unit Jacobian, rate identity, inverse transport: {roots} roots')
    # Independent symbolic check of logarithmic kernel exponents, no floating point.
    t, rate, parent_rate, d = s.symbols('t rate parent_rate d', positive=True)
    zero(s.limit((rate/eps-1)*s.log(1-eps*t), eps, 0, dir='+') + rate*t)
    zero(s.limit(-(parent_rate/eps)*s.log(1+eps*d*(1-eps*t)/t), eps, 0, dir='+')
         + parent_rate*d/t)
    print('Exact logarithmic kernel limits: 2')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    check()
