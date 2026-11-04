#!/usr/bin/env python3
"""Exact algebraic checks for Tkachev, arXiv:2610.13532v1, Problem 11.3(b).

This program checks polynomial certificates, not the analytic directional
limit or an unbounded family by enumeration. The accompanying preregistration
supplies the general argument. No floating point is used.
"""
from itertools import combinations_with_replacement
from math import factorial, prod
import sys
import sympy as S


def zero(expr):
    assert S.simplify(S.expand(expr)) == 0, expr


def laplacian(p, xs):
    return S.expand(sum(S.diff(p, x, 2) for x in xs))


def infinity_laplacian(p, xs):
    g = S.Matrix([S.diff(p, x) for x in xs])
    return S.expand((g.T * S.hessian(p, xs) * g)[0])


def fischer(p, q, xs):
    pp, qq = S.Poly(p, *xs), S.Poly(q, *xs)
    return S.expand(sum(c * qq.coeff_monomial(e) * prod(factorial(a) for a in e)
                        for e, c in pp.terms()))


def sphere_integral(p, xs):
    """Normalized sphere monomial moments, in dimension len(xs)."""
    out = S.Integer(0)
    for exps, c in S.Poly(S.expand(p), *xs).terms():
        if any(e % 2 for e in exps):
            continue
        half_degree = sum(exps) // 2
        numerator = prod(prod(range(1, e, 2)) for e in exps)
        denominator = prod(len(xs) + 2*j for j in range(half_degree))
        out += c * S.Rational(numerator, denominator)
    return S.simplify(out)


def block(xs, m):
    return xs[0] * (sum(x*x for x in xs[1:1+m]) - sum(x*x for x in xs[1+m:]))


def positive_nodes(m, blocks):
    d, n = 2*m+1, blocks*(2*m+1)
    alpha, beta = S.sqrt(3)/3, S.sqrt(6)/3
    nodes = []
    for b in range(blocks):
        for j in range(m):
            for sign in [-1, 1]:
                for negative in [False, True]:
                    v = [S.Integer(0)] * n
                    v[b*d] = -alpha if negative else alpha
                    v[b*d+1+(m if negative else 0)+j] = sign * beta
                    nodes.append(tuple(v))
    return nodes


def check_sup_certificates():
    A, B, Y, Z = S.symbols('A B Y Z', nonnegative=True)
    zero(4*(A+B)**3 - 27*A*B**2 - (2*A-B)**2*(A+4*B))
    zero(A*(Y+Z)**2 - A*(Y-Z)**2 - 4*A*Y*Z)
    print('SUP_CERTIFICATES exact factorization and nonnegative difference OK')


def check_fidelity():
    xs = S.symbols('v0 y1 y2 y3 y4')
    u = block(xs, 2)/2
    assert laplacian(u, xs) == 0
    zero(infinity_laplacian(u, xs) - 2*sum(x*x for x in xs)*u)
    ff = fischer(u, u, xs)
    l2 = sphere_integral(u*u, xs)
    assert ff == 2 and l2 == S.Rational(2, 315)
    nodes = positive_nodes(2, 1)
    for a in nodes:
        assert sum(c*c for c in a) == 1
        zero(u.subs(dict(zip(xs, a))) - S.sqrt(3)/9)
    ratio = S.Rational(1, 27) / l2
    assert ratio == S.Rational(35, 6)
    # Lemma 6.1, including its explicit x1*x2*x3 normalization test.
    h = xs[0]*xs[1]*xs[2]
    assert fischer(h, h, xs) == 1
    assert sphere_integral(h*h, xs) == S.Rational(1, 315)
    print(f'FIDELITY Pol(1,4), source labels pol14 and tab:ratio: Fischer={ff}, L2sq={l2}, sup2=1/27, R={ratio}')
    print('FIDELITY harmonic Hsiang identity: infinity_laplacian(u)=2*|x|^2*u OK')


def check_family(m):
    d, n = 2*m+1, 4*m+2
    xs = S.symbols(f'x0:{n}')
    q1, q2 = block(xs[:d], m), block(xs[d:], m)
    p = q1 + q2
    r1, r2 = sum(x*x for x in xs[:d]), sum(x*x for x in xs[d:])
    assert S.Poly(p, *xs).total_degree() == 3
    assert all(sum(e) == 3 for e, _ in S.Poly(p, *xs).terms())
    assert laplacian(p, xs) == 0
    infp = infinity_laplacian(p, xs)
    zero(infp - 8*(r1*q1+r2*q2))
    nodes = positive_nodes(m, 2)
    assert len(nodes) == 8*m and len(set(nodes)) == 8*m
    M = 2*S.sqrt(3)/9
    for a in nodes:
        zero(sum(c*c for c in a) - 1)
        zero(p.subs(dict(zip(xs, a))) - M)
    moment = S.expand(sum(sum(x*c for x, c in zip(xs, a))**3 for a in nodes) / len(nodes))
    zero(moment - p/(2*m*S.sqrt(3)))
    ff = fischer(p, p, xs)
    l2 = sphere_integral(p*p, xs)
    N = n*(n+2)*(n+4)
    assert ff == 8*m and l2 == S.Rational(8*m, N)
    # Check the whole cubic dual space, including mixed-block perturbations.
    checked = 0
    for indices in combinations_with_replacement(range(n), 3):
        h = prod(xs[i] for i in indices)
        average = sum(prod(a[i] for i in indices) for a in nodes) / len(nodes)
        fh = fischer(p, h, xs)
        zero(average - fh/(12*m*S.sqrt(3)))
        zero(sphere_integral(p*h, xs) - fh/N)
        checked += 1
    assert checked == n*(n+1)*(n+2)//6
    # A rationally parameterized infinite subset of the maximum set.
    t = S.symbols('t', real=True)
    xi1, xi2 = (1-t*t)/(1+t*t), 2*t/(1+t*t)
    zero(xi1*xi1+xi2*xi2-1)
    zero(xi2/(1+xi1)-t)
    a = [S.Integer(0)]*n
    a[0], a[1], a[2] = S.sqrt(3)/3, S.sqrt(6)/3*xi1, S.sqrt(6)/3*xi2
    zero(sum(c*c for c in a)-1)
    zero(p.subs(dict(zip(xs, a)))-M)
    # Two exact evaluations exclude every real lambda: use X if lambda != 8,
    # and Y if lambda == 8. The original Hsiang theta is -lambda.
    X, Y = [S.Integer(0)]*n, [S.Integer(0)]*n
    X[0] = X[1] = Y[0] = Y[1] = Y[d] = S.Integer(1)
    readings = []
    lam = S.symbols('lambda', real=True)
    residuals = []
    for point in [X, Y]:
        sub = dict(zip(xs, point))
        pv, normsq, iv = p.subs(sub), sum(c*c for c in point), infp.subs(sub)
        readings.append((pv, normsq, iv))
        residuals.append(S.expand(iv-lam*normsq*pv))
    assert readings == [(1, 2, 16), (1, 3, 16)]
    assert residuals == [16-2*lam, 16-3*lam]
    assert S.solve(residuals, [lam]) == []
    # Kernel of v -> D_v p is zero, so every ambient coordinate is essential.
    gradients = [S.Poly(S.diff(p, x), *xs) for x in xs]
    monomials = sorted(set(e for g in gradients for e, _ in g.terms()))
    matrix = S.Matrix([[g.coeff_monomial(e) for g in gradients] for e in monomials])
    assert matrix.rank() == n
    print(f'FAMILY m={m}, n={n}: homogeneous, harmonic, full infinity-Laplacian identity OK')
    print(f'CUBATURE nodes={len(nodes)}, mixed-direction cubic monomials={checked}: complete third moment and Fischer/sphere identities OK')
    print(f'NORMS Fischer={ff}, L2sq={l2}, sup2=4/27, R={S.Rational(4,27)/l2}')
    print(f'EXTREMAL infinite rational parameterization OK; annihilating-direction rank={n}')
    print(f'NON_HSIANG exact (p,|x|^2,Delta_inf) readings={readings}; incompatible lambda=8 and 16/3 OK')


def main():
    print(f'Python {sys.version.split()[0]}; SymPy {S.__version__}; exact algebra only')
    check_sup_certificates()
    check_fidelity()
    for m in [2, 3, 4]:
        check_family(m)
    print('SCOPE: the directional-derivative argument and general-m proof are prose obligations, not numerical conclusions.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
