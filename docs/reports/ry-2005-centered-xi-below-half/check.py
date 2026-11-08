#!/usr/bin/env python3
"""Exact checks of the centred-xi obstruction; numerical diagnostics are labelled.

Requires Python 3, SymPy and mpmath. No external program source is used.
The exact checks establish algebraic steps, not the analytic zero-sum theorem.
"""
from fractions import Fraction
import sys
import sympy as s
import mpmath as mp


def exact_checks():
    # Source fidelity: (V.21) and ξ(s)=ξ(1-s), at the regular pair -1, 2.
    def xi_at(x):
        return s.simplify(x * (x - 1) * s.pi ** (-x / 2)
                          * s.gamma(x / 2) * s.zeta(x) / 2)
    assert xi_at(s.Integer(2)) == s.pi / 6
    assert xi_at(s.Integer(-1)) == s.pi / 6
    print('EXACT fidelity: xi(2) = xi(-1) = pi/6 (V.21 and reflection)')

    a, y = s.symbols('a y', real=True)
    numerator = a**4 - 6*a*a*y*y + y**4
    inverse_fourth_real = s.re(s.expand_complex((a-s.I*y)**(-4)))
    assert s.factor(inverse_fourth_real - numerator/(a*a+y*y)**4) == 0
    print('EXACT conjugate-zero kernel: Re((a-iy)^(-4)) = (a^4-6a^2*y^2+y^4)/(a^2+y^2)^4')

    # Continuous family certificate. w=1-a²>=0 and v=y²-196>0.
    # Every coefficient of P(1-w,196+v) is a nonnegative integer.
    w, v = s.symbols('w v', nonnegative=True)
    certificate = s.Poly(s.expand((1-w)**2 - 6*(1-w)*(196+v)
                                 + (196+v)**2), w, v)
    assert all(c >= 0 for c in certificate.coeffs())
    assert certificate.coeff_monomial(1) == 37241
    print('EXACT family certificate:', certificate.as_expr(), '>= 37241 > 0')
    for n in range(2, 33):
        for k in range(1, n):
            alpha = Fraction(k, 2*n)
            assert 0 < alpha < Fraction(1,2)
            for beta in [Fraction(0), Fraction(1,2), Fraction(1)]:
                aa = alpha-beta
                yy = Fraction(14) + Fraction(1,n)
                assert aa**4-6*aa**2*yy**2+yy**4 > 0
    print('EXACT rational spot checks: every k/(2n), 2<=n<=32, 1<=k<n; beta=0,1/2,1')

    u = s.symbols('u', real=True)
    A, C = s.symbols('A C', real=True)
    B = s.symbols('B', positive=True)
    f = s.exp(-A*u*u - s.sqrt(2)*B*u**3/3 - C*u**4/6)
    jet = s.series(f, u, 0, 5).removeO().expand()
    assert jet.coeff(u,2) == -A
    assert jet.coeff(u,3) == -s.sqrt(2)*B/3
    difference = s.expand(jet.subs(u,s.sqrt(2)*u)-2*jet+1)
    coeff = s.simplify(difference.coeff(u,3))
    assert s.simplify(coeff + (4-2*s.sqrt(2))*B/3) == 0
    assert coeff.is_negative is True
    assert all(s.simplify(difference.coeff(u,j)) == 0 for j in range(3))
    print('EXACT small-t second difference: -(4-2*sqrt(2))*g3/3 * t^(3/2) + O(t^2), g3>0')

    z = s.symbols('z', real=True)
    assert s.expand(z*z-2*z+1-(1-z)**2) == 0
    print('EXACT Laplace obstruction: z^2-2z+1=(1-z)^2>=0 for real z')


def numerical_diagnostics():
    # These evaluations are diagnostics only; they are not exact evidence.
    mp.mp.dps = 100
    def xi(x):
        if x == 0 or x == 1:
            return mp.mpf(1)/2
        return x*(x-1)*mp.power(mp.pi,-x/2)*mp.gamma(x/2)*mp.zeta(x)/2
    alpha = mp.mpf(1)/4
    g = lambda x: mp.log(xi(x))
    derivatives = [mp.diff(g,alpha,n) for n in range(1,5)]
    g1,g2,g3,g4 = derivatives
    print('NUMERICAL ONLY: mpmath dps=100, alpha=1/4')
    for n,value in enumerate(derivatives,1):
        print('NUMERICAL g'+str(n)+' = '+mp.nstr(value,45))
    assert g3 > 0 and g4 < 0
    def F(t):
        lam = mp.sqrt(2*t)
        return mp.exp(g(alpha)-g(alpha+lam)+g1*lam)
    asymptotic = -(4-2*mp.sqrt(2))*g3/3
    print('NUMERICAL predicted limit = '+mp.nstr(asymptotic,45))
    for power in [4,6,8,10,12,16,20]:
        t = mp.mpf(10)**(-power)
        delta = F(2*t)-2*F(t)+1
        scaled = delta/t**mp.mpf('1.5')
        assert delta < 0
        print('NUMERICAL t=10^(-'+str(power)+'): second_difference='
              +mp.nstr(delta,35)+'; scaled='+mp.nstr(scaled,35))
    assert abs(scaled/asymptotic-1) < mp.mpf('0.000001')


if __name__ == '__main__':
    print('Python '+sys.version.split()[0]+'; sympy '+s.__version__+'; mpmath '+mp.__version__)
    exact_checks()
    numerical_diagnostics()
    print('ALL_CHECKS_OK')
