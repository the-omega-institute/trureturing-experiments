#!/usr/bin/env python3
"""Exact fidelity/coefficient checks; separate, explicitly non-certified diagnostics.

Source: Hake--Keller--Pogorzelski, arXiv:2608.24213v1, section 4,
Remark after Proposition 10. No Lean proof is attempted.
"""
import argparse
import functools
import itertools
import sympy as S


def exact_checks():
    d, z = S.symbols('d z')
    p = 1 - d / 2
    a = d * (d - 2) * ((d + 2) * z - 3) / 48
    d4 = p * (p - 2) * (3*d + 6*(p-4) + (p-4)*(p-6)*z)
    derived = 2*(d+1)*a/d - d4/(24*d)
    source = (d-2)**2/(16*d) * (5*(d+2)**2*z/8 - (2*d+3))
    assert S.factor(derived-source) == 0
    print('Q_difference =', S.factor(derived-source))
    lap = (p-2)*(p+d-4) - 4*(d+2)
    assert S.factor(lap + (d-2)**2/4 + 4*(d+1)) == 0
    print('laplacian_coefficient_difference = 0')
    # Check the inverse-symbol coefficients without assigning a dimension.
    t, rho, A4, A6, A8 = S.symbols('t rho A4 A6 A8', nonzero=True)
    denom = (t**2*rho/2 - t**4*A4/24 + t**6*A6/720 - t**8*A8/40320)/d
    kept = 2*d/(t**2*rho) + d*A4/(6*rho**2) + t**2*(d*A4**2/(72*rho**3)-d*A6/(180*rho**2))
    assert all(S.expand((denom*kept).series(t,0,6).removeO()).coeff(t,k) == (1 if k == 0 else 0) for k in range(6))
    print('reciprocal_symbol_coefficients = OK')
    # Fidelity: source Lemma 11, normalized covariance identity, for SRW.
    # This is symbolic in all coordinates, not a sample of x-values.
    for n in range(3,9):
        x = S.symbols('x:'+str(n)); rr = sum(v*v for v in x)
        lhs = sum(S.Rational(1,2*n)*(sign*x[j])**2/rr for j in range(n) for sign in (-1,1))
        assert S.factor(lhs-S.Rational(1,n)) == 0
        h = sum(v**4 for v in x)-3*rr**2/S.Integer(n+2)
        assert S.expand(sum(S.diff(h,v,2) for v in x)) == 0
        # Independently differentiate a radial power at nonzero x for fixed d.
        k = 1-S.Rational(n,2)
        real_d4 = sum(S.diff(rr**(k/2),v,4) for v in x)
        expected = k*(k-2)*rr**((k-4)/2)*(3*n+6*(k-4)+(k-4)*(k-6)*sum(v**4 for v in x)/rr**2)
        assert S.simplify(real_d4-expected) == 0
    print('FIDELITY: Lemma 11 SRW covariance identity = 1/d, d=3,...,8, symbolic coordinates')
    print('harmonic_quartic_and_direct_fourth_derivatives = OK, d=3,...,8')
    # Rational direction tests, covering every nonnegative small vector up to symmetry.
    count = 0
    for n in range(3,7):
        for x in itertools.combinations_with_replacement(range(4),n):
            rr = sum(v*v for v in x)
            if not rr: continue
            zz = S.Rational(sum(v**4 for v in x),rr**2)
            assert S.factor((derived-source).subs({d:n,z:zz})) == 0
            count += 1
    print('small_lattice_direction_exact_cases =',count)
    assert S.factor(source.subs(z,1)-(d-2)**2*(5*d**2+4*d-4)/(128*d)) == 0
    assert S.factor(source.subs(z,1/d)+(d-2)**2*(11*d**2+4*d-20)/(128*d**2)) == 0
    print('axis_and_diagonal_sign_formulas = OK')


def numerical_checks():
    # These floating-point checks are diagnostics, not exact evidence or a proof
    # of the uniform remainder. Independent Green values use the heat integral
    # G(x)= integral_0^infinity exp(-t) product_j I_|x_j|(t/d) dt.
    import mpmath as mp
    mp.mp.dps = 42

    def scaled_bessel(n,z):
        if z < max(10000,100*(n+1)**2):
            return mp.exp(-z)*mp.besseli(n,z)
        term = mp.mpf(1); total = term
        for k in range(1,41):
            term *= -mp.mpf(4*n*n-(2*k-1)**2)/(8*k*z)
            total += term
        return total/mp.sqrt(2*mp.pi*z)

    @functools.lru_cache(None)
    def green(x):
        n=len(x); rr=sum(v*v for v in x)
        def f(u):
            if not u:
                return 2*n*rr/(2*mp.pi*rr)**(mp.mpf(n)/2) if n == 3 else mp.mpf(0)
            zz=mp.mpf(rr)/(u*u)
            return 2*n*rr/u**3 * mp.fprod(scaled_bessel(v,zz) for v in x)
        return mp.quad(f,[0,1,2,4,8,16,mp.inf],method='tanh-sinh',maxdegree=9)

    def G(x): return green(tuple(sorted(abs(v) for v in x)))
    def residual(x):
        n=len(x); rr=sum(v*v for v in x); g=G(x)
        mean=mp.mpf(0); green_mean=mp.mpf(0)
        for j in range(n):
            for sign in (-1,1):
                y=list(x);y[j]+=sign
                gy=G(tuple(y))
                mean += mp.sqrt(gy/g)/(2*n)
                green_mean += gy/(2*n)
        harmonic_error=abs(green_mean/g-1)
        assert harmonic_error < mp.mpf("1e-25"), (x,harmonic_error)
        w=1-mean
        q=mp.mpf(n-2)**2/(16*n)*(mp.mpf(5)*(n+2)**2/8*mp.mpf(sum(v**4 for v in x))/rr**2-(2*n+3))
        return (w-mp.mpf(n-2)**2/(8*n*rr)-q/rr**2)*rr**3

    print('NUMERICAL_DIAGNOSTICS_ONLY: 42 decimal digits; 40-term Bessel tail; no certified error bound',flush=True)
    for n in (3,4):
        for direction in [(1,)+(0,)*(n-1),(1,1)+(0,)*(n-2),(1,)*n]:
            values=[]
            for scale in (4,8,12):
                x=tuple(scale*v for v in direction); val=residual(x)
                assert mp.isfinite(val) and abs(val) < 100
                values.append(mp.nstr(val,16))
                print('d='+str(n),'direction='+str(direction),'scale='+str(scale),'r^6_residual='+mp.nstr(val,16),flush=True)
            print('sample_max_abs < 100; finite samples do not establish an asymptotic bound',flush=True)

if __name__ == '__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--exact-only',action='store_true');args=parser.parse_args()
    exact_checks()
    if not args.exact_only: numerical_checks()
    print('ALL_CHECKS_OK')
