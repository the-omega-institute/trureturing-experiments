#!/usr/bin/env python3
"""Exact checks for arXiv:2411.10439v2, Conjecture 5.1 (Python 3 + SymPy).

Checks finite identities and proposed saddle algebra. Does not certify the
analytic limit, uniform errors or tails. Uses no floating point or quadrature.
"""
from math import factorial
import sympy as s

Q = s.Rational


def equal(a, b, label):
    difference = s.simplify(s.expand_log(a - b, force=True))
    if difference != 0:
        difference = difference.replace(
            lambda e: e.is_Pow, lambda e: s.Pow(s.factor(e.base), e.exp))
        difference = s.simplify(difference)
    assert difference == 0, (label, difference)
    assert not (s.sympify(a).atoms(s.Float) | s.sympify(b).atoms(s.Float)), label


def source_fidelity():
    # Source label KintAsympt1: integrate the half-integer K expansion
    # using integral_0^infinity x^j exp(-x) dx = j!.
    count = 0
    for n in range(1, 20, 2):
        q = (n + 1) // 2  # nu = n/2 + 1 = q + 1/2
        for m in range(11):
            exact_integral = s.sqrt(s.pi / 2) * sum(
                Q(factorial(q+j), factorial(j)*factorial(q-j)*2**j)
                * factorial(2*m+q-j) for j in range(q+1))
            source = (s.sqrt(s.pi)*s.gamma(Q(n+3, 2))*2**Q(n, 2)
                      * Q(factorial(m+q)*factorial(2*m), factorial(m)*factorial(q)))
            equal(exact_integral, source, ('source K integral', n, m))
            count += 1
    print('FIDELITY: source KintAsympt1, odd n=1..19, m=0..10:', count, 'exact equalities')

    # Source label convexQuadratic: independent Gaussian integral and powers.
    y = s.symbols('y', real=True)
    count = 0
    for p in (Q(1, 4), Q(1, 2), s.Integer(1), s.Integer(2), s.Integer(4)):
        integral = s.integrate(s.exp(-p*y*y/(2*(1+p))), (y, -s.oo, s.oo))
        base = 2*s.pi*s.sqrt((1+p)**((1+p)/p)/p)
        computed = s.sqrt(2*s.pi)*(1+p)**(1/(2*p))*integral
        for n in range(1, 9):
            equal(computed**n, base**n, ('source Gaussian Mahler', p, n))
            count += 1
    print('FIDELITY: source convexQuadratic, five exact p, n=1..8:', count, 'exact equalities')


def finite_ball_identities():
    # Exact functional identities for all r>0 in the small odd dimensions.
    # This compares the whole integral with Bessel I, not just series terms.
    t, r = s.symbols('t r', positive=True)
    for n in range(1, 8, 2):
        a = s.gamma(1+Q(n, 2))/(s.sqrt(s.pi)*s.gamma(Q(n+1, 2)))
        marginal = a*s.integrate(s.exp(r*t)*(1-t*t)**((n-1)//2), (t, -1, 1))
        source = s.gamma(1+Q(n, 2))*(2/r)**Q(n, 2)*s.besseli(Q(n, 2), r)
        difference = s.cancel(s.expand((s.simplify(marginal)-s.simplify(source)).rewrite(s.exp)))
        assert difference == 0, ('whole ball marginal', n, difference)
        assert not (marginal.atoms(s.Float) | source.atoms(s.Float))
    print('FINITE IDENTITIES: whole ball marginal equals normalized Bessel I for all r>0, n=1,3,5,7: 4')
    # Beta-integral marginal moments versus coefficients in the source's
    # Gamma(1+n/2) (2/r)^(n/2) I_(n/2)(r), label hpB2n.
    count = 0
    for n in range(1, 13):
        alpha = Q(n-1, 2)
        normalization = s.gamma(Q(n+2, 2))/(s.sqrt(s.pi)*s.gamma(Q(n+1, 2)))
        for k in range(13):
            moment = normalization*s.gamma(Q(2*k+1, 2))*s.gamma(alpha+1)/s.gamma(k+alpha+Q(3, 2))
            marginal_coefficient = moment/factorial(2*k)
            bessel_coefficient = s.gamma(1+Q(n, 2))/(4**k*factorial(k)*s.gamma(k+1+Q(n, 2)))
            equal(marginal_coefficient, bessel_coefficient, ('ball coefficient', n, k))
            if n % 2 == 1:
                t = s.symbols('t')
                direct = normalization*s.integrate(t**(2*k)*(1-t*t)**int(alpha), (t, -1, 1))
                equal(direct, moment, ('direct polynomial moment', n, k))
            count += 1
    print('FINITE IDENTITIES: ball marginal versus source Bessel series, n=1..12, k=0..12:', count)
    print('FINITE IDENTITIES: odd-dimensional moments also checked by polynomial integration: 78')

    # Corollary 3.8: exact change r=n*z in the radial normalization.
    count = 0
    for n in range(1, 13):
        v = s.pi**Q(n, 2)/s.gamma(1+Q(n, 2))
        for p in (Q(1, 4), Q(1, 2), s.Integer(1), s.Integer(2), s.Integer(4)):
            radial = n*v*v*(s.Integer(n)/p)**n
            transformed = radial/s.Integer(n)**n/(s.gamma(1+Q(n, 2))*2**Q(n, 2))**(1/p)
            source = n*s.pi**n/(2**(Q(n, 2)/p)*p**n*s.gamma(1+Q(n, 2))**(2+1/p))
            equal(transformed, source, ('radial Corollary 3.8', n, p))
            count += 1
    print('FINITE IDENTITIES: radial normalization versus Corollary 3.8:', count)


def proposed_asymptotic_algebra():
    t = s.symbols('t', positive=True)
    z = t/(1-t*t)
    h = z*t + s.log(1-t*t)/2
    hprime = s.diff(h, t)/s.diff(z, t)
    hsecond = s.diff(hprime, t)/s.diff(z, t)
    equal(hprime, t, 'inner envelope derivative')
    equal(hsecond, (1-t*t)**2/(1+t*t), 'inner envelope curvature')
    p = s.symbols('p', positive=True)
    tsq = p/(1+p)
    zs = s.sqrt(p*(1+p))
    hs = p - s.log(1+p)/2
    fs = s.log(p*(1+p))/2 - 1 + s.log(1+p)/(2*p)
    equal(zs*s.sqrt(tsq), p, 'outer critical point')
    equal(s.log(zs)-hs/p, fs, 'outer phase value')
    curvature = 1/zs**2 + ((1-tsq)**2/(1+tsq))/p
    equal(curvature, 2/(p*(1+2*p)), 'outer curvature')
    amplitude = s.sqrt((1-tsq)/(1+tsq))
    equal(amplitude, 1/s.sqrt(1+2*p), 'inner amplitude')
    base = 2*s.pi*s.sqrt((1+p)**((1+p)/p)/p)
    equal(s.log(2*s.pi*s.E/p)+fs, s.log(base), 'exponential rate')
    computed = (1/zs)*amplitude**(-1/p)*s.sqrt(2*s.pi/curvature)
    proposed = s.sqrt(s.pi)*(1+2*p)**((1+p)/(2*p))/s.sqrt(1+p)
    equal(computed, proposed, 'outer prefactor')
    # n*v_n^2 contributes 1/pi to the non-exponential prefactor.
    equal(proposed/s.pi, (1+2*p)**((1+p)/(2*p))/s.sqrt(s.pi*(1+p)), 'final prefactor')
    print('EXPECTED SETTLEMENT: saddle, exponential base and prefactor identities hold symbolically for p>0')
    print('SCOPE: no convergence, uniform Laplace error or tail estimate is certified by this program')


def main():
    print('Exact arithmetic only; SymPy', s.__version__)
    source_fidelity()
    finite_ball_identities()
    proposed_asymptotic_algebra()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
