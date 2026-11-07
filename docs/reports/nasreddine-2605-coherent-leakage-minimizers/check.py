#!/usr/bin/env python3
"""Exact fidelity and rigidity checks for arXiv:2605.11024v2.

These are finite algebraic checks and symbolic identities, not a proof of the
quantified entropy minimizer classification. Requires Python 3 and SymPy.
"""
from itertools import product
from fractions import Fraction
import sympy as s


def xlog(x):
    return s.S.Zero if x == 0 else x * s.log(x)


def phi(a, e, x):
    radical = s.sqrt((a - e) ** 2 + 4 * x)
    lp, lm = (a + e + radical) / 2, (a + e - radical) / 2
    return xlog(lp) + xlog(lm) - xlog(a) - xlog(e)


def fidelity():
    # The source's Asymptotic sharpness proposition, pure two-level family: spectrum {0,1}, entropy h(q).
    count = 0
    for den in range(3, 13):
        for num in range(1, den):
            if 2 * num >= den:
                continue
            q = s.Rational(num, den)
            m = s.Matrix([[1 - q, s.sqrt(q * (1 - q))],
                          [s.sqrt(q * (1 - q)), q]])
            assert m.trace() == 1 and m.det() == 0 and m * m == m
            assert s.simplify(phi(1 - q, q, q * (1 - q)) +
                              xlog(q) + xlog(1 - q)) == 0
            count += 1
    print(f"FIDELITY pure-family spectrum and exact entropy identity: {count} cases")


def symbolic_core():
    t, a, e, x, delta = s.symbols('t a e x delta', positive=True)
    lp, lm = s.symbols('lambda_plus lambda_minus', positive=True)
    # The integrand uses only trace and determinant of the active block.
    assert s.expand((t + a) * (t + e) - x -
                    (t*t + (a+e)*t + a*e-x)) == 0
    primitive = sum(sign * ((t+u)*s.log(t+u)-(t+u))
                    for sign, u in [(1,lp),(1,lm),(-1,a),(-1,e)])
    assert s.expand(s.diff(primitive,t) -
                    (s.log(t+lp)+s.log(t+lm)-s.log(t+a)-s.log(t+e))) == 0
    monotone_gap = x/((a+t)*(e+t)) - x/((a+delta+t)*(e+t))
    assert s.factor(monotone_gap - x*delta/((a+t)*(a+delta+t)*(e+t))) == 0
    e1, e2, z1, z2 = s.symbols('e1 e2 z1 z2')
    gap = (t+z1+z2)*(t+e1)*(t+e2)-(t+z1)*(t+z2)*(t+e1+e2)
    decomposition = t*(e1*e2-z1*z2)+e1*z1*(e2-z2)+e2*z2*(e1-z1)
    assert s.expand(gap-decomposition) == 0
    u, v = s.symbols('u v', nonnegative=True)
    # z_i=e_i-loss_i; first summand already forces strictness when
    # both leakages are positive and one loss is positive.
    assert s.expand(e1*e2-(e1-u)*(e2-v)-(u*e2+(e1-u)*v)) == 0
    print('SYMBOLIC_CORE determinant, primitive derivative, monotonicity and merging identities: OK')


def scalar_grid():
    grid = [Fraction(i,4) for i in range(5)]
    times = [Fraction(1,3), Fraction(1), Fraction(2)]
    pairs = 0
    for e1, e2 in product(grid, repeat=2):
        for f1, f2 in product(grid, repeat=2):
            z1, z2 = e1*(1-f1), e2*(1-f2)
            for t in times:
                left = (t+z1)*(t+z2)/((t+e1)*(t+e2))
                right = (t+z1+z2)/(t+e1+e2)
                assert left <= right
                strict = e1 > 0 and e2 > 0 and (z1 < e1 or z2 < e2)
                assert (left < right) == strict
                pairs += 1
    # Complete enumeration of this explicitly bounded scalar data set.
    # Each tuple is (a_j,e_j,x_j); no floating-point entropy comparison.
    abar = Fraction(1)
    entries = [(a,e,a*e*f) for a in [abar/2,abar]
               for e in [Fraction(0),Fraction(1,4),Fraction(1,2)]
               for f in [Fraction(0),Fraction(1,2),Fraction(1)]]
    collapse, equalities = 0, 0
    for n in range(1,4):
        for data in product(entries, repeat=n):
            total_x = sum(x for a,e,x in data)
            if not total_x:
                continue
            for rem in [Fraction(0),Fraction(1,4)]:
                total_e = rem + sum(e for a,e,x in data)
                positive = [j for j,(_,e,_) in enumerate(data) if e > 0]
                expected_equal = (rem == 0 and len(positive) == 1 and
                                  data[positive[0]][0] == abar)
                for t in times:
                    original = Fraction(1)
                    for a,e,x in data:
                        original *= 1-x/((a+t)*(e+t))
                    merged = 1-total_x/((abar+t)*(total_e+t))
                    assert 0 < original <= merged
                    assert (original == merged) == expected_equal
                    collapse += 1
                    equalities += original == merged
    print(f'SCALAR_GRID two-term merging: {pairs} exact comparisons')
    print(f'SCALAR_GRID full collapse, n=1..3: {collapse} comparisons; {equalities} equalities, all single-channel')


def unit_vector(dim, rotated):
    v = s.zeros(dim,1)
    if dim == 1 or not rotated:
        v[0] = 1
    else:
        v[0],v[1] = s.Rational(3,5),s.Rational(4,5)
    assert (v.conjugate().T*v)[0] == 1
    return v


def displayed_family():
    z = s.symbols('z')
    count = singular = floor_saturated = 0
    for dp,dq in product(range(1,4),repeat=2):
        a0 = s.Rational(1,2*dp+2)
        for e in [s.Rational(1,4),s.Rational(3,8),1-dp*a0]:
            astar = 1-e-(dp-1)*a0
            assert astar >= a0 and e > 0
            for theta in [s.Rational(1,4),s.Rational(1,2),s.S.One]:
                c = theta*astar*e
                for phase in [s.S.One,s.I,(3+4*s.I)/5]:
                    assert s.expand(phase*s.conjugate(phase)) == 1
                    p,q = unit_vector(dp,True),unit_vector(dq,True)
                    pp,qq = p*p.T,q*q.T
                    a = astar*pp+a0*(s.eye(dp)-pp)
                    b = s.sqrt(c)*phase*p*q.T
                    cc = e*qq
                    rho = a.row_join(b).col_join(b.conjugate().T.row_join(cc))
                    assert rho == rho.conjugate().T and rho.trace() == 1
                    assert cc.trace() == e
                    assert s.simplify((b.conjugate().T*b).trace()-c) == 0
                    # Characteristic polynomials certify all eigenvalues exactly:
                    # active roots have nonnegative trace and determinant;
                    # spectators are a0 or zero. This is a PSD certificate.
                    expected = (z-a0)**(dp-1)*z**(dq-1)*(z*z-(astar+e)*z+astar*e-c)
                    assert s.expand(rho.charpoly(z).as_expr()-expected) == 0
                    assert astar*e-c >= 0 and astar+e > 0
                    assert s.expand((a-a0*s.eye(dp)).charpoly(z).as_expr()-
                                    z**(dp-1)*(z-(astar-a0))) == 0
                    # Pinching spectra: same protected spectators; leakage
                    # spectators zero. Cancelling them gives exactly Phi.
                    pinch = s.diag(a,cc)
                    expected_pinch = (z-a0)**(dp-1)*z**(dq-1)*(z-astar)*(z-e)
                    assert s.expand(pinch.charpoly(z).as_expr()-expected_pinch) == 0
                    count += 1
                    singular += theta == 1
                    floor_saturated += astar == a0
    print(f'DISPLAYED_FAMILY equation (14), dp,dq=1..3: {count} exact certificates; {singular} singular endpoints; {floor_saturated} saturated floors')


def main():
    fidelity()
    symbolic_core()
    scalar_grid()
    displayed_family()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
