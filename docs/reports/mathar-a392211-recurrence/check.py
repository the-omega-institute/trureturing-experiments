#!/usr/bin/env python3
"""Exact certificates and source-fidelity checks for Mathar's A392211 recurrence."""
import hashlib
import json
import re
from pathlib import Path

import sympy as sp

ROOT = Path(__file__).resolve().parent
x, y, n, t = sp.symbols('x y n t')


def zero(expression):
    assert sp.Poly(expression, x, y, domain=sp.QQ).is_zero


def read_table():
    result = {}
    for line in (ROOT / 'b392211.txt').read_text().splitlines():
        if line.strip() and not line.lstrip().startswith('#'):
            index, value = map(int, line.split())
            assert index not in result
            result[index] = value
    assert sorted(result) == list(range(1, 1031))
    return result


def source_polynomials(source):
    # Parse the actual conjecture text, independently of the supplied certificate.
    line = next(s for s in source['formula'] if s.startswith('Conjecture D-finite'))
    expr = line.split('recurrence ', 1)[1].split('=0.', 1)[0]
    expr = re.sub(r'a\(n(?:-(\d+))?\)',
                  lambda m: f'z{m.group(1) or 0}', expr)
    z = sp.symbols('z0:6')
    expr = sp.sympify(expr.replace('^', '**'), locals=dict(zip(map(str, z), z)))
    p = [sp.expand(expr).coeff(zj) for zj in z]
    assert sp.expand(expr - sum(pj*zj for pj, zj in zip(p, z))) == 0
    return p


def check_certificate(source, cert, p):
    P = 25*y**4 - 35*y**3 + (7+20*x)*y**2 + (3-14*x)*y - 3*x + 4*x**2
    zero(P - sp.sympify(cert['P']))
    for actual, supplied in zip(p, cert['p_j'], strict=True):
        assert sp.expand(actual - sp.sympify(supplied)) == 0
    Q, T = sp.diff(P, y), sp.diff(P, x)
    assert Q.subs({x: 0, y: 0}) == 3
    N, powers = [y, -x*T], [0, 1, 3, 5, 7]
    for r in powers[1:4]:
        last = N[-1]
        N.append(sp.expand(x*((sp.diff(last, x)*Q-r*last*sp.diff(Q, x))*Q -
                              (sp.diff(last, y)*Q-r*last*sp.diff(Q, y))*T)))
    for actual, supplied in zip(N, cert['N_k'], strict=True):
        zero(actual-sp.sympify(supplied))
    assert powers == cert['denominator_powers']
    op = sp.Poly(sum(x**j*pj.subs(n, t+j) for j, pj in enumerate(p)), t)
    assert op.degree() == 4
    c = [op.nth(k) for k in range(5)]
    for actual, supplied in zip(c, cert['c_k'], strict=True):
        zero(actual-sp.sympify(supplied))
    forcing = -151822515600*x**4 + 171245697000*x**3
    zero(forcing-sp.sympify(cert['forcing']))
    zero(forcing-sp.Rational(7, 20)*c[0])
    H = sp.Poly(sp.sympify(cert['polynomial_quotient']), y, x, domain=sp.ZZ)
    lhs = sp.Poly(sum(c[k]*N[k]*Q**(7-powers[k]) for k in range(5)) -
                  forcing*Q**7, y, x, domain=sp.ZZ)
    quotient, remainder = lhs.div(sp.Poly(P, y, x, domain=sp.ZZ))
    assert remainder.is_zero and quotient == H
    zero(lhs.as_expr()-P*H.as_expr())
    assert len(H.terms()) == 136
    print('EXACT_SYZYGY_OK: integer polynomial remainder 0; H has 136 terms')
    zero(c[4]+209*(1140304*x-572535)*(6400*x**4-9419*x**2+256))
    assert c[4].subs(x, 0) == 30632912640
    theta = lambda f: x*sp.diff(f, x)
    zero(theta(theta(forcing))-7*theta(forcing)+12*forcing)
    print('DFINITE_CERTIFICATE_OK: (theta-3)(theta-4) kills the forcing; c4 is nonzero')

    # Algebraic bridge to FORMULA (3.a) and the square-root branch in NAME/(3.b).
    B = (4*x-7*y+10*y**2)/3
    zero(B**2-B+y**2-y-sp.Rational(4, 9)*P)
    back = (-4*x-7*B+10*B**2)/3
    zero(back-y-sp.Rational(40, 27)*P)
    _, rem = sp.div(sp.Poly(P.subs({x: -x, y: B}, simultaneous=True), y, x),
                    sp.Poly(P, y, x))
    assert rem.is_zero
    print('SOURCE_BRANCH_ALGEBRA_OK: reflection and square-root identities modulo P')


def check_series(source, table, p):
    bound = max(table)
    a, square, cube, fourth = ([0]*(bound+1) for _ in range(4))
    # Triangular quartic coefficient recursion, independent of the conjectured recurrence.
    for k in range(1, bound+1):
        square[k] = sum(a[j]*a[k-j] for j in range(1, k))
        cube[k] = sum(a[j]*square[k-j] for j in range(1, k))
        fourth[k] = sum(square[j]*square[k-j] for j in range(1, k))
        rhs = (3*(k == 1)-4*(k == 2)+14*a[k-1]-20*square[k-1] -
               7*square[k]+35*cube[k]-25*fourth[k])
        value, remainder = divmod(rhs, 3)
        assert remainder == 0
        a[k] = value
        assert a[k] == table[k], ('source mismatch', k)
    initial = list(map(int, source['data'].split(',')))
    assert a[1:len(initial)+1] == initial
    assert a[1:13] == [1, 1, 5, 11, 80, 207, 1710, 4755, 41900, 121364, 1112045, 3306877]
    for k in range(bound+1):
        sign = (-1)**k
        # Source FORMULA (3.a), literally after coefficient extraction.
        assert 3*a[k] == -4*(k == 1)-7*sign*a[k]+10*sign*square[k]
        # S=1-2A is the square root of 1+4A(-x)-4A(-x)^2; S(0)=1.
        assert -4*a[k]+4*square[k] == 4*sign*a[k]-4*sign*square[k]
        assert 25*fourth[k]-35*cube[k]+7*square[k]+20*(square[k-1] if k else 0) +\
            3*a[k]-14*(a[k-1] if k else 0)-3*(k == 1)+4*(k == 2) == 0
    print(f'SOURCE_FIDELITY_OK: quartic, NAME square-root branch, FORMULA (3.a), and {bound} b-file terms')
    integer_polys = [list(map(int, sp.Poly(pj, n).all_coeffs())) for pj in p]
    def evaluate(coefficients, k):
        value = 0
        for coefficient in coefficients:
            value = value*k+coefficient
        return value
    def residual(k):
        return sum(evaluate(pj, k)*a[k-j] for j, pj in enumerate(integer_polys))
    for k in range(6, bound+1):
        assert residual(k) == 0, ('recurrence mismatch', k)
        assert evaluate(integer_polys[0], k) > 0
    assert a[0] == 0 and residual(5) == 0
    # Positivity for every integer n>=6, by a polynomial with positive shifted coefficients.
    leading_shift = sp.Poly(p[0].subs(n, t+6), t)
    assert all(value > 0 for value in leading_shift.all_coeffs())
    print(f'RECURRENCE_OK: {bound-5} checks, n=6..{bound}; n=5 also holds with a(0)=0')
    print('LEADING_COEFFICIENT_OK: p0(n)>0 for every integer n>=6')
    print('FIRST_TERMS:', ','.join(map(str, a[1:16])))


def main():
    source = json.loads((ROOT / 'source.json').read_text())
    assert source['number'] == 392211 and source['revision'] == 35 and source['offset'] == '1,3'
    cert = json.loads((ROOT / 'certificate.json').read_text())
    p = source_polynomials(source)
    check_certificate(source, cert, p)
    check_series(source, read_table(), p)
    print('BFILE_SHA256:', hashlib.sha256((ROOT / 'b392211.txt').read_bytes()).hexdigest())
    print('ARITHMETIC: integers and SymPy exact rational polynomials; no floating point')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
