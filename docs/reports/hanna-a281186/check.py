#!/usr/bin/env python3
"""Exact checks for Hanna's OEIS A281186 real-exponent positivity conjecture.

Python standard library only. No floating point or network access is used.
The finite sign computation is evidence for a later Lean refutation, not a
kernel-checked theorem. The interval conclusion also needs the analytic proof
and Pringsheim's theorem in the accompanying preregistration.
"""
from fractions import Fraction as F
from hashlib import sha256
from math import comb, factorial

PUBLISHED = [
    1, 3, 15, 132, 1485, 22248, 389259, 8252280, 195922881,
    5374450440, 161160732999, 5426241074712, 196597702744245,
    7846633439714088, 333286781355259875, 15383602537188087672,
    749475751265692473321, 39275782157858062883400,
    2158712482065435156322719, 126614778653174063332276440,
]
EXPECTED_HASH = 'effab943471a3e14c2c2c740a3375b3c5bf4195151ce57e51529f4d8f70f8acd'


def mul(a, b):
    """Ordinary power-series product, truncated to the common list length."""
    assert len(a) == len(b)
    return [sum((a[k] * b[n-k] for k in range(n+1)), F(0))
            for n in range(len(a))]


def inv(a):
    assert a[0] != 0
    b = [1 / a[0]]
    for n in range(1, len(a)):
        b.append(-sum((a[k]*b[n-k] for k in range(1, n+1)), F(0))/a[0])
    return b


def exp(a):
    assert a[0] == 0
    b = [F(1)]
    for n in range(1, len(a)):
        b.append(sum((k*a[k]*b[n-k] for k in range(1, n+1)), F(0))/n)
    return b


def log(a):
    assert a[0] == 1
    reciprocal = inv(a)
    # The n-th coefficient of log(a) is [x^(n-1)] a'/a divided by n.
    return [F(0)] + [sum((k*a[k]*reciprocal[n-k]
                            for k in range(1, n+1)), F(0))/n
                     for n in range(1, len(a))]


def primitive(a):
    return [F(0)] + [a[n-1]/n for n in range(1, len(a))]


def source_step(g, t):
    power = exp([t*c for c in log(g)])
    return exp(mul(inv(g), primitive(power)))


def direct_source_solution(t, degree):
    # Integration gains one coefficient, so degree iterations determine the jet.
    g = [F(1)] + [F(0)]*degree
    for _ in range(degree):
        g = source_step(g, t)
    assert source_step(g, t) == g
    return g


def numerator_rows(degree):
    # Q_n = (((5n+3)z+6-n) Q_(n-1)) + 3(1+z) Q'_(n-1).
    # All entries are signed Python integers, with no totalized subtraction.
    row = [1]
    constants = [1]
    for n in range(1, degree+1):
        old = row
        row = []
        for k in range(n+1):
            a = old[k] if k < len(old) else 0
            b = old[k-1] if k > 0 else 0
            c = old[k+1] if k+1 < len(old) else 0
            row.append((6-n+3*k)*a + (5*n+3)*b + 3*(k+1)*c)
        constants.append(row[0])
    return constants


def lagrange_check(N):
    # z=y/3, w=x/3, w=H(z), H'=(1+3z)exp(-5z), h=H/z.
    # h_j = j! [z^j]h. With F=exp(3z)h^(-N), use h F'=(3h-Nh')F.
    h = [F(1)] + [F((-5)**(j-1)*(3*j-5), j+1) for j in range(1, N)]
    f = [F(1)]
    for n in range(N-1):
        rhs = sum((comb(n,k)*(3*h[k]-N*h[k+1])*f[n-k]
                   for k in range(n+1)), F(0))
        lhs_known = sum((comb(n,k)*h[k]*f[n-k+1]
                         for k in range(1,n+1)), F(0))
        f.append(rhs-lhs_known)
    assert f[-1].denominator == 1
    return f[-1].numerator


def eval_poly(p, z):
    return sum((c*z**k for k,c in enumerate(p)), F(0))


def check_derivative_numerators():
    # Independently differentiate the symbolic rational shape, at exact rational
    # z values. This checks signs and index shifts in the chain-rule numerator.
    p = [1]
    for n in range(1, 13):
        q = [0]*(len(p)+1)
        for k,c in enumerate(p):
            q[k] += (6-n)*c
            q[k+1] += (5*n+3)*c
            if k:
                q[k-1] += 3*k*c
                q[k] += 3*k*c
        dp = [k*p[k] for k in range(1,len(p))]
        for z in [F(-3,4), F(0), F(1,5), F(2)]:
            expected = ((5*n+3)*(1+z)-3*(2*n-1))*eval_poly(p,z)
            expected += 3*(1+z)*eval_poly(dp,z)
            assert eval_poly(q,z) == expected
        p = q


def main():
    assert F(8,3) < F(267,100)
    g = direct_source_solution(F(8,3), 19)
    terms = [g[n]*factorial(n)*3**n for n in range(20)]
    assert terms == PUBLISHED
    print('Fidelity: literal exp/inverse/zero-primitive iteration matches all 20 source terms.')
    print('Published a(0)..a(19):', ','.join(map(str, PUBLISHED)))
    # Independently check the differential equation stated in FORMULA, to degree 18.
    power = exp([F(8,3)*c for c in log(g)])
    deriv = [(n+1)*g[n+1] for n in range(19)] + [F(0)]
    one_log = log(g)
    one_log[0] += 1
    assert mul(deriv,one_log)[:19] == power[:19]
    print('Fidelity: A^(8/3) = A\'*(1+log A), coefficients 0..18, exact.')
    check_derivative_numerators()
    constants = numerator_rows(719)
    assert [1] + [3*constants[n-1] for n in range(1,20)] == PUBLISHED
    q = constants[719]
    assert q < 0
    assert all(c > 0 for c in constants[:719])
    signed = str(q)
    digest = sha256(signed.encode('ascii')).hexdigest()
    assert digest == EXPECTED_HASH
    independent = lagrange_check(720)
    assert independent == q
    print('Witness: Q_719(0) =', q)
    print('Signed-integer ASCII SHA-256:', digest)
    print('Independent Lagrange recurrence: f_719 = Q_719(0), exact equality.')
    print('Sign: [x^720]G_(8/3) = Q_719(0)/(3^719*720!) < 0.')
    print('Earlier degrees: [x^m]G_(8/3) > 0 for 0 <= m < 720 (checked).')
    # This rational majorant bounds exp(167/100) for the entire real interval.
    a = F(167,100)
    partial = sum((a**j/factorial(j) for j in range(13)), F(0))
    tail_majorant = (a**13/factorial(13))/(1-a/14)
    M = partial+tail_majorant
    assert M == F(40786381559745290820669898287650613487,
                  7677916646400000000000000000000000000)
    assert a/14 < 1
    assert M < F(5313,1000) < F(16,3)
    print('Interval Taylor majorant M =', M)
    print('Exact inequalities: M < 5313/1000 < 16/3 = 2*(8/3).')
    print('Family inputs: exp(a)>1+a for a>0 and monotonicity, together with this')
    print('majorant, give 0 < |x_minus(t)| < R_plus(t) on [8/3,267/100].')
    print('Infinitely many negatives on that interval requires the analytic inverse')
    print('and Pringsheim argument; the arithmetic checks do not certify it alone.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
