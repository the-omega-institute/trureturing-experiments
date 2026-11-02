#!/usr/bin/env python3
"""Exact source-fidelity and separated-star counterexample certificates.

Original code; no floating-point arithmetic. Requires Python 3 and SymPy.
Ideal intersections use elimination; ideal quotients use (J intersect (z))/z.
"""
from itertools import combinations
import sys
import sympy as sp


def gb(generators, variables, domain):
    return sp.groebner(generators, *variables, order="lex", domain=domain)


def member(basis, polynomial):
    return basis.reduce(sp.expand(polynomial))[1] == 0


def same_ideal(left, right, variables, domain):
    a, b = gb(left, variables, domain), gb(right, variables, domain)
    assert all(member(a, p.as_expr()) for p in b.polys)
    assert all(member(b, p.as_expr()) for p in a.polys)


def intersection(left, right, variables, domain):
    t = sp.Dummy("elimination_parameter")
    basis = gb([t * p for p in left] + [(1 - t) * p for p in right],
               (t,) + tuple(variables), domain)
    return [p.as_expr() for p in basis.polys if not p.as_expr().has(t)]


def colon(left, z, variables, domain):
    answer = []
    for p in intersection(left, [z], variables, domain):
        q, rem = sp.div(p, z, *variables, domain=domain)
        assert rem == 0
        answer.append(q)
    return answer


def minimal_covers(vertices, edges):
    covers = []
    for size in range(len(vertices) + 1):
        for subset in combinations(vertices, size):
            c = frozenset(subset)
            if all(c.intersection(e) for e in edges) and not any(a <= c for a in covers):
                covers.append(c)
    return set(covers)


def eval_zero(p, assignment, variables, domain):
    return sp.Poly(sp.expand(p.subs(assignment, simultaneous=True)),
                   *variables, domain=domain).is_zero


def fidelity(domain):
    # Source Example 3.6, using Example 3.1's six-vertex graph.
    a, b, c, d, e, f = variables = sp.symbols("a b c d e f")
    edges = [(a, b), (b, c), (b, e), (d, e), (e, f)]
    ideal = [u * v for u, v in edges]
    q1, q2 = frozenset([b, e]), frozenset([a, c, e])
    covers = minimal_covers(variables, edges)
    assert q1 in covers and q2 in covers
    j = ideal + [a + b]
    basis = gb(j, variables, domain)
    assert not member(basis, b)
    same_ideal(colon(j, b, variables, domain), [a, b, c, e], variables, domain)
    # The minimal-prime calculation after a=-b is monomial and exact.
    reduced_vertices = (b, c, d, e, f)
    reduced = [-b * b, b * c, b * e, d * e, e * f]
    assert all(eval_zero(p, {a: -b}, reduced_vertices, domain)
               for p in [u.subs(a, -b) - v for u, v in zip(ideal, reduced)])
    reduced_covers = minimal_covers(reduced_vertices,
                                    [(b,), (b, c), (b, e), (d, e), (e, f)])
    assert frozenset([b, e]) in reduced_covers
    # (a,b,e) is prime, properly contained in the associated prime (a,b,c,e).
    assert not member(gb([a, b, e], variables, domain), c)
    print(f"FIDELITY {domain}: Example 3.6, J:b=(a,b,c,e); (a,b,e) minimal; embedded extension verified")


def family_case(r, s, domain):
    xs = sp.symbols(f"x1:{r + 1}")
    ys = sp.symbols(f"y1:{s + 1}")
    variables = xs + ys
    edges = [(x, y) for x in xs for y in ys]
    ideal = [x * y for x, y in edges]
    linear = sum(ys)
    star = xs[0] + linear
    j = ideal + [star]
    basis = gb(j, variables, domain)
    # Definition 2.1: every edge containing the center contains a listed leaf;
    # every listed leaf has maximum generator exponent exactly one.
    assert all(xs[0] not in edge or any(y in edge for y in ys) for edge in edges)
    assert all(max(int(y in edge) for edge in edges) == 1 for y in ys)
    assert len({xs[0], *ys}) == s + 1
    # Singleton collections are pairwise separated. The graph is connected,
    # and each vertex has a neighbor; no isolated-vertex convention is used.
    assert all(any(v in edge for edge in edges) for v in variables)
    assert minimal_covers(variables, edges) == {frozenset(xs), frozenset(ys)}
    same_ideal(intersection(list(xs), list(ys), variables, domain), ideal, variables, domain)
    # Exact epsilon-map certificate in k[e]/(e^2), including characteristic two.
    epsilon = sp.Symbol("epsilon")
    assignment = {v: 0 for v in variables}
    assignment[xs[0]], assignment[ys[0]] = epsilon, -epsilon
    eps_basis = gb([epsilon**2], (epsilon,), domain)
    assert all(member(eps_basis, p.subs(assignment, simultaneous=True)) for p in j)
    assert not member(eps_basis, xs[0].subs(assignment))
    # Independent Groebner and ideal-quotient checks of the socle certificate.
    assert not member(basis, xs[0])
    assert all(member(basis, v * xs[0]) for v in variables)
    same_ideal(colon(j, xs[0], variables, domain), list(variables), variables, domain)
    # Coefficient-one polynomial identities explain the characteristic-free mechanism.
    assert all(sp.expand(x * xs[0] - (x * star - sum(x * y for y in ys))) == 0 for x in xs)
    # Neither extension is m: maps to the base field kill its generators
    # and send an omitted variable to one.
    map_a = {v: 0 for v in variables}
    map_a[ys[0]], map_a[ys[1]] = 1, -1
    assert all(eval_zero(p, map_a, variables, domain) for p in list(xs) + [star])
    assert sp.Poly(ys[0].subs(map_a), *variables, domain=domain).as_expr() == 1
    map_b = {v: 0 for v in variables}
    map_b[xs[1]] = 1
    assert all(eval_zero(p, map_b, variables, domain) for p in list(ys) + [star])
    assert sp.Poly(xs[1].subs(map_b), *variables, domain=domain).as_expr() == 1
    # f is regular modulo I: exact colon equality, independently of the socle check.
    same_ideal(colon(ideal, star, variables, domain), ideal, variables, domain)
    # Structural cross-check after eliminating x1=-L. Conca-Herzog Lemma 3.2
    # gives the general primary decomposition; here its equality and its
    # three irredundancy witnesses are checked by exact elimination.
    reduced_variables = xs[1:] + ys
    product = [x * y for x in xs[1:] for y in ys] + [linear * y for y in ys]
    n_squared = [u * v for i, u in enumerate(reduced_variables) for v in reduced_variables[i:]]
    t_ideal = list(xs[1:]) + [linear]
    y_inter_t = intersection(list(ys), t_ideal, reduced_variables, domain)
    decomposition = intersection(y_inter_t, n_squared, reduced_variables, domain)
    same_ideal(product, decomposition, reduced_variables, domain)
    by, bt, bn = (gb(ig, reduced_variables, domain) for ig in (list(ys), t_ideal, n_squared))
    assert member(by, linear) and member(bt, linear) and not member(bn, linear)
    assert member(bt, xs[1]**2) and member(bn, xs[1]**2) and not member(by, xs[1]**2)
    assert member(by, ys[0]**2) and member(bn, ys[0]**2) and not member(bt, ys[0]**2)
    print(f"FAMILY r={r} s={s} {domain}: star; J:x1=m; x1!=0; old primes=A,B; m!=A+(f),B+(f); I:f=I; YT=Y intersect T intersect n^2")


def main():
    print(f"Python {sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}; SymPy {sp.__version__}")
    domains = [sp.QQ, sp.GF(2), sp.GF(3), sp.GF(5)]
    for domain in domains:
        fidelity(domain)
        for r in range(2, 5):
            for s in range(2, 5):
                family_case(r, s, domain)
    print("SCOPE: 4 source-fidelity cases; 36 family instances (2<=r,s<=4), exact Q/F2/F3/F5 arithmetic")
    print("BOUNDARY: finite checks support the supplied characteristic-free family argument; they are not a Lean proof")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
