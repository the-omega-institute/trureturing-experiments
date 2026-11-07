#!/usr/bin/env python3
"""Exact finite audits of Bouthat--Chavez--Sheng Section 5(ii), arXiv:2606.17356v1.

These computations test the proposed structural answer; they are not a Lean proof.
All mathematical evidence uses integers, Fraction, or exact SymPy algebra.
"""
from fractions import Fraction as Q
from itertools import combinations, combinations_with_replacement
from math import comb, factorial
import sympy as sp


def moment(law, k):
    return sum((w * x**k for x, w in law), Q(0))


def centered_data(law):
    assert sum(w for _, w in law) == 1
    mean = moment(law, 1)
    centered = [(x - mean, w) for x, w in law]
    var = moment(centered, 2)
    assert var > 0
    t = moment(centered, 4) / var**2
    a = moment(centered, 3)**2 / var**3
    u = moment(centered, 6) / var**3
    assert t >= 1 and 0 <= a <= t - 1 and u >= t*t + a
    return centered, var, t, a, u


def direct_moment(coeffs, law, d):
    # Independently convolve the actual iid atomic sums, without cumulants.
    masses = {Q(0): Q(1)}
    for c in coeffs:
        new = {}
        for x, w in masses.items():
            for y, v in law:
                z = x + c*y
                new[z] = new.get(z, Q(0)) + w*v
        masses = new
    return sum((w*x**d for x, w in masses.items()), Q(0))


def power_sums(coeffs):
    return tuple(sum(x**j for x in coeffs) for j in (2, 3, 4, 6))


def coefficients(ps):
    p2, p3, p4, p6 = ps
    return p6, 15*(p2*p4-p6), 10*(p3*p3-p6), 15*p2**3-45*p2*p4+30*p6


def F(ps, t, a, u):
    A, B, D, C = coefficients(ps)
    return A*u+B*t+D*a+C


def S(ps, t):
    return F(ps, t, Q(0), t*t)


def T(ps, t):
    return F(ps, t, t-1, t*t+t-1)


def laws():
    result = []
    for size, patterns in [(2, [(1, 1), (1, 2)]), (3, [(1, 1, 1), (1, 1, 2), (1, 2, 3)])]:
        for support in combinations(range(-2, 3), size):
            for weights in patterns:
                total = sum(weights)
                result.append([(Q(x), Q(w, total)) for x, w in zip(support, weights)])
    result.append([(Q(x), Q(1, 5)) for x in (-2, -1, 0, 1, 3)])
    return result


def trace_powers(n, mask):
    edges = list(combinations(range(n), 2))
    A = [[0]*n for _ in range(n)]
    for i, (v, w) in enumerate(edges):
        if mask & (1 << i):
            A[v][w] = A[w][v] = 1
    B = [[int(i == j) for j in range(n)] for i in range(n)]
    ps = {}
    for k in range(1, 7):
        B = [[sum(B[i][r]*A[r][j] for r in range(n)) for j in range(n)] for i in range(n)]
        ps[k] = sum(B[i][i] for i in range(n))
    assert ps[1] == 0
    return tuple(ps[k] for k in (2, 3, 4, 6))


def symbolic_checks():
    t, a, u, s, x, z = sp.symbols('t a u s x z', real=True)
    A, B, D, C = sp.symbols('A B D C', real=True)
    theta = a/(t-1)
    f = A*u+B*t+D*a+C
    sym = A*t*t+B*t+C
    two = A*t*t+(A+B+D)*t+C-A-D
    assert sp.cancel(f-((1-theta)*sym+theta*two+(u-t*t-a)*A)) == 0
    # Verify the two squared-polynomial expectation identities symbolically.
    m = {0: 1, 1: 0, 2: 1, 3: s, 4: t, 5: sp.Symbol('m5'), 6: u}
    def expectation(poly):
        return sp.expand(sum(coef*m[monomial[0]] for monomial, coef in sp.Poly(poly, x).terms()))
    assert sp.expand(expectation((x*x-s*x-1)**2)-(t-s*s-1)) == 0
    assert sp.expand(expectation((x**3-t*x-s)**2)-(u-t*t-s*s)) == 0
    # Actual two-point law: verify probability normalization and moments through six.
    root = sp.sqrt(s*s+4)
    rp, rm = (s+root)/2, (s-root)/2
    wp, wm = -rm/root, rp/root
    assert sp.simplify(wp+wm) == 1
    for k, expected in [(1, 0), (2, 1), (3, s), (4, s*s+1), (6, s**4+3*s*s+1)]:
        assert sp.simplify(wp*rp**k+wm*rm**k-expected) == 0
    # Actual symmetric law: weights 1/(2t), 1-1/t, 1/(2t).
    tp = sp.symbols('tp', positive=True)
    for k, expected in [(1, 0), (2, 1), (3, 0), (4, tp), (6, tp**2)]:
        assert sp.simplify(((-sp.sqrt(tp))**k+sp.sqrt(tp)**k)/(2*tp)-expected) == 0
    AA, BB, CC, DD, EE, FF = sp.symbols('AA BB CC DD EE FF')
    P, R = AA*t*t+BB*t+CC, DD*t*t+EE*t+FF
    derivative = (AA*EE-BB*DD)*t*t+2*(AA*FF-CC*DD)*t+(BB*FF-CC*EE)
    assert sp.expand(sp.diff(P, t)*R-P*sp.diff(R, t)-derivative) == 0
    print('SYMBOLIC: square identities, boundary moments, conical decomposition, derivative OK')


def fidelity():
    rad = [(Q(-1), Q(1, 2)), (Q(1), Q(1, 2))]
    count = 0
    for n in range(2, 9):
        for d in (2, 4, 6, 8):
            complete = (n-1,)+(-1,)*(n-1)
            matching = (1, -1)*(n//2)+(0,)*(n % 2)
            # Proposition in Section 3.5, label prop:SpecialNormsRademacher, (a),(b).
            expected_k = Q(2)**(d+1-n)*sum(comb(n-1, k)*k**d for k in range(n))/factorial(d)
            expected_m = Q(2)**d*sum(comb(2*(n//2), k)*(n//2-k)**d for k in range(2*(n//2)+1))/(Q(4)**(n//2)*factorial(d))
            assert direct_moment(complete, rad, d)/factorial(d) == expected_k
            assert direct_moment(matching, rad, d)/factorial(d) == expected_m
            count += 2
        ps = power_sums(complete)
        p2, _, p4, p6 = ps
        assert direct_moment(complete, rad, 6) == 15*p2**3-30*p2*p4+16*p6
    for r in range(1, 7):
        for s in range(1, 7):
            for d in (2, 4, 6, 8):
                # K_{r,s} has nonzero eigenvalues +/-sqrt(rs): d is even.
                actual = (r*s)**(d//2)*direct_moment((1, -1), rad, d)/factorial(d)
                assert actual == Q((4*r*s)**(d//2), 2*factorial(d))
                count += 1
    print(f'FIDELITY: {count} exact source Rademacher comparator evaluations; sixth-power formula OK')


def expansion_audit():
    all_laws = laws()
    count = positive_residual = zero_residual = 0
    lists = []
    for n in range(2, 6):
        for coeffs in combinations_with_replacement(range(-2, 3), n):
            if sum(coeffs) == 0 and any(coeffs):
                lists.append(coeffs)
                ps = power_sums(coeffs)
                for law in all_laws:
                    centered, var, t, a, u = centered_data(law)
                    actual = direct_moment(coeffs, law, 6)/var**3
                    assert actual == F(ps, t, a, u)
                    assert actual == direct_moment(coeffs, centered, 6)/var**3
                    count += 1
                    positive_residual += int(u > t*t+a)
                    zero_residual += int(u == t*t+a)
    assert positive_residual and zero_residual
    print(f'IID EXPANSION: {len(lists)} trace-zero integer multisets of length 2..5, entries -2..2; {len(all_laws)} laws; {count} exact checks')
    print(f'RESIDUAL: positive={positive_residual}, zero={zero_residual}')
    return all_laws


def graph_audit(all_laws):
    data = [centered_data(law)[2:] for law in all_laws]
    graph_count = pair_count = 0
    for n in range(2, 6):
        complete = power_sums((n-1,)+(-1,)*(n-1))
        matching = power_sums((1, -1)*(n//2)+(0,)*(n % 2))
        comparators = [complete, matching, (2, 0, 2, 2), (12, 24, 84, 732)]
        for mask in range(1, 1 << (n*(n-1)//2)):
            G = trace_powers(n, mask)
            graph_count += 1
            for t, a, u in data:
                fg = F(G, t, a, u)
                assert fg >= G[0]**3 > 0
                if t == 1:
                    assert a == 0 and u == 1
                for H in comparators:
                    fh = F(H, t, a, u)
                    sh, th = S(H, t), T(H, t)
                    assert fh >= H[0]**3 > 0 and sh > 0 and th > 0
                    ratio = fg/fh
                    if t > 1:
                        theta, residual = a/(t-1), u-t*t-a
                        w = ((1-theta)*sh/fh, theta*th/fh, residual*H[3]/fh)
                        r = (S(G, t)/sh, T(G, t)/th, Q(G[3], H[3]))
                        assert all(v >= 0 for v in w) and sum(w) == 1
                        assert ratio == sum(wi*ri for wi, ri in zip(w, r))
                        assert min(r) <= ratio <= max(r)
                    else:
                        assert ratio == S(G, 1)/S(H, 1)
                    pair_count += 1
    print(f'GRAPHS: all {graph_count} labeled nonempty-edge graphs of order 2..5; {pair_count} exact graph/comparator/law checks')


def interior_optimizer():
    t = sp.symbols('t', real=True)
    G, H = (6, 0, 18, 54), (12, 24, 84, 732)
    assert trace_powers(4, 7) == G  # Edge order 01,02,03,12,13,23.
    assert trace_powers(4, 63) == H
    sg, sh = S(G, t), S(H, t)
    tg, th = T(G, t), T(H, t)
    ds = sp.factor(sp.diff(sg, t)*sh-sg*sp.diff(sh, t))
    dt = sp.factor(sp.diff(tg, t)*th-tg*sp.diff(th, t))
    assert sp.expand(ds-19440*(-19*t*t+14*t+105)) == 0
    assert sp.expand(dt+58320*(t+3)**2) == 0
    star = (7+2*sp.sqrt(511))/19
    other = (7-2*sp.sqrt(511))/19
    assert sp.simplify(star-1).is_positive and other.is_negative
    value = sp.simplify((sg/sh).subs(t, star))
    rad, limit = sp.Rational(9, 77), sp.Rational(9, 122)
    assert sp.simplify((sg/sh).subs(t, 1)) == rad
    assert sp.limit(sg/sh, t, sp.oo) == limit
    assert sp.limit(tg/th, t, sp.oo) == limit
    assert sp.simplify(value-rad).is_positive and sp.simplify(value-limit).is_positive
    # Establish derivative signs exactly; no decimal approximations.
    assert sp.expand(ds+19440*(19*t*t-14*t-105)) == 0
    assert sp.simplify(ds.subs(t, star)) == 0
    print('INTERIOR: G=K_{1,3}, H=K_4; tuples', G, H)
    print('S derivative numerator =', ds)
    print('T derivative numerator =', dt)
    print('unique global maximum at t* =', star)
    print('maximum sixth-power ratio =', value, '> 9/77 > 9/122')


if __name__ == '__main__':
    symbolic_checks()
    fidelity()
    all_laws = expansion_audit()
    graph_audit(all_laws)
    interior_optimizer()
    print('ALL_CHECKS_OK')
