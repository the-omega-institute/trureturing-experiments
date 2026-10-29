#!/usr/bin/env python3
"""Exact certificates for the weighted kernel in Acosta Babb, 2410.17681v1.

Only integers and Fraction enter the mathematical checks. Analytic bridges
(alternating-series theorem, continuity, Sturm comparison, integral formula,
and Taylor's theorem) are explained in README.md and are not Lean proofs.
"""
from fractions import Fraction as Q
from math import factorial
import json

N = 40
A = (Q(24048255, 10**7), Q(24048256, 10**7))
B = (Q(55200781, 10**7), Q(55200782, 10**7))
U = (Q(5191, 10000), Q(5192, 10000))
V = (Q(3402, 10000), Q(3404, 10000))


def term(nu, k, z):
    return (z / 2) ** (2*k + nu) / (factorial(k)*factorial(k+nu))


def enclosure(nu, z):
    """Alternating bracket P_40 <= J_nu <= P_40 + t_40 for z >= 0.

    For k >= 40, t_(k+1)/t_k = (z/2)^2/((k+1)(k+nu+1)) <= 1.
    Terms tend to zero (indeed the ratios are eventually <= 1/2).
    """
    assert z >= 0 and (z/2)**2 < (N+1)*(N+nu+1)
    partial = sum(((-1)**k * term(nu, k, z) for k in range(N)), Q(0))
    tail = term(nu, N, z)
    assert tail < Q(1, 10**55)
    return partial, partial + tail


def scaled(q):
    """Outward interval on a 10^-12 grid, for readable exact output."""
    den = 10**12
    return [str(Q((q[0]*den).__floor__(), den)),
            str(Q((q[1]*den).__ceil__(), den))]


def add(i, j):
    return (i[0]+j[0], i[1]+j[1])


def mul(i, j):
    products = [x*y for x in i for y in j]
    return (min(products), max(products))


def polynomial_interval(nu, z):
    """Termwise outward interval plus uniform alternating-tail bracket."""
    out = (Q(0), Q(0))
    for k in range(N):
        lo, hi = term(nu, k, z[0]), term(nu, k, z[1])
        out = add(out, (lo, hi) if k % 2 == 0 else (-hi, -lo))
    return (out[0], out[1] + term(nu, N, z[1]))


def sqrt_interval(z):
    den = 10**12
    # Integer binary search; no float, Decimal, or approximate algebra.
    low, high = 0, 10*den
    while high-low > 1:
        mid = (low+high)//2
        if Q(mid, den)**2 <= z:
            low = mid
        else:
            high = mid
    assert Q(low, den)**2 <= z < Q(high, den)**2
    return Q(low, den), Q(high, den)


def main():
    # Fidelity: Eq. (2.8) in the journal / eqn:intJ0rdr in the TeX.
    # Integrating each term of y*J0(y) gives precisely the coefficient
    # of z*J1(z); J0'=-J1 is checked coefficient by coefficient too.
    for k in range(N):
        integrated = Q((-1)**k, 2**(2*k)*factorial(k)**2*(2*k+2))
        z_j1 = Q((-1)**k, 2**(2*k+1)*factorial(k)*factorial(k+1))
        assert integrated == z_j1
        derivative = Q((-1)**(k+1)*(2*k+2), 2**(2*k+2)*factorial(k+1)**2)
        assert derivative == -z_j1
    print('FIDELITY: integral_0^z y*J0(y)dy = z*J1(z), 40 exact coefficient identities')
    print('FIDELITY: J0_prime = -J1, 40 exact coefficient identities')

    signs = {}
    for name, interval, sign in [('a', A, 1), ('b', B, -1)]:
        lo, hi = map(lambda z: enclosure(0, z), interval)
        assert (lo[0] > 0 and hi[1] < 0) if sign == 1 else (lo[1] < 0 and hi[0] > 0)
        signs[name] = [scaled(lo), scaled(hi)]
    print('ROOT_ENDPOINT_SIGNS ' + json.dumps(signs, sort_keys=True))
    # The first two root indices follow from positivity on [0,2] and
    # Sturm spacing on [2,6]: 4*pi/sqrt(17)>3, using pi>31/10.
    assert Q(16)*Q(31, 10)**2 > 9*17
    assert A[0] > 2 and B[1] < 6 and B[1]-A[0] < 6
    print('ROOT_ISOLATION_INPUTS: 2<a<3, 5<b<6; spacing>3 excludes extra roots')

    # Transfer J1 values from the lower rational endpoint to the root:
    # |J1'|<=1; the uncertainty is at most the root-interval width.
    ea, eb = enclosure(1, A[0]), enclosure(1, B[0])
    aw, bw = A[1]-A[0], B[1]-B[0]
    root_u = (ea[0]-aw, ea[1]+aw)
    root_v = (-eb[1]-bw, -eb[0]+bw)
    assert U[0] < root_u[0] < root_u[1] < U[1]
    assert V[0] < root_v[0] < root_v[1] < V[1]
    print('J1_ROOT_BOUNDS ' + json.dumps({'J1(a)': scaled(root_u), '-J1(b)': scaled(root_v)}, sort_keys=True))
    assert A[1]**7 < Q(342, 100)**5
    assert B[0]**7 > Q(1093, 100)**5
    derivative_lower = (Q(1093,100)*V[0] - 2*Q(342,100)*U[1])/3
    assert derivative_lower == Q(27843,500000) > Q(1,20)
    print('UNIFORM_DERIVATIVE_LOWER ' + str(derivative_lower))
    # Uniform for alpha in [2/5,1/2], since a,b>1 and b/a>1.
    assert A[1] < Q(2405,1000) and B[1] < Q(5521,1000)
    assert Q(1551,1000)**2 > A[1]
    assert Q(235,100)**2 > B[1]
    curvature_upper = (2*Q(2405,1000)**2*Q(1551,1000) + Q(5521,1000)**2*Q(235,100))/6
    assert curvature_upper < 15
    hmax = Q(1,1000)
    assert -Q(1,20) + Q(15,2)*hmax < 0
    print('UNIFORM_CURVATURE_UPPER ' + str(curvature_upper))
    print('UNIFORM_FAMILY: K_2^alpha(1-h) <= h*(-1/20 + (15/2)*h) < 0 for 2/5<=alpha<=1/2, 0<h<=1/1000')

    # alpha_* is defined by (b/a)^(alpha_*+1)*v=2*u.
    # Rational exponent comparisons become exact integer powers.
    assert B[1]**67 * V[1]**50 < A[0]**67 * (2*U[0])**50
    assert B[0]**27 * V[0]**20 > A[1]**27 * (2*U[1])**20
    print('SHARP_THRESHOLD: 17/50 < alpha_star < 7/20')

    # Direct enclosure of the two J0 values at a specified x, independent
    # of the endpoint derivative/Taylor route.  This is a genuine x<1.
    x = Q(999,1000)
    ja = polynomial_interval(0, (A[0]*x, A[1]*x))
    jb = polynomial_interval(0, (B[0]*x, B[1]*x))
    assert ja[0] > 0 and jb[1] < 0
    sa = (sqrt_interval(A[0])[0], sqrt_interval(A[1])[1])
    sb = (sqrt_interval(B[0])[0], sqrt_interval(B[1])[1])
    k = mul((Q(1,3),Q(1,3)), add(mul((Q(2),Q(2)),mul(sa,ja)),mul(sb,jb)))
    assert k[1] < 0
    print('DIRECT_POINT ' + json.dumps({'alpha':'1/2','M':2,'x':str(x),'J0(a*x)':scaled(ja),'J0(b*x)':scaled(jb),'K':scaled(k)},sort_keys=True))
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
