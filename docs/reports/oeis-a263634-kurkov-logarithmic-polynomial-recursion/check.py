#!/usr/bin/env python3
"""Exact finite checks of OEIS A263634 revision 166 and Kurkov's recursion.

Requires Python 3 and SymPy. These checks do not prove the unbounded claim.
The target ring is QQ[x1,...,x8]; no variable is numerically specialized
for the polynomial-identity or quotient-family checks.
"""
from functools import lru_cache
from math import comb, factorial
import sympy as sp

N = 8
xs = sp.symbols('x1:9')
ZERO = sp.S.Zero
ONE = sp.S.One


def multiply_truncated(a, b):
    return [sp.expand(sum((a[j] * b[d-j] for j in range(d+1)), ZERO))
            for d in range(N+1)]


def logarithmic_polynomials():
    # Coefficients of U = sum x_i*t^i/i!, and the finite log(1+U).
    u = [ZERO] + [xs[i-1] * sp.Rational(1, factorial(i))
                  for i in range(1, N+1)]
    power = [ONE] + [ZERO] * N
    log_coeffs = [ZERO] * (N+1)
    for r in range(1, N+1):
        power = multiply_truncated(power, u)
        scalar = sp.Rational((-1)**(r+1), r)
        log_coeffs = [sp.expand(a + scalar*b)
                      for a, b in zip(log_coeffs, power)]
    return [sp.expand(factorial(i)*log_coeffs[i]) for i in range(N+1)]


@lru_cache(None)
def R(n, k):
    assert n >= 1 and k >= 1 and n+k-1 <= N
    if n == 1:
        return xs[k-1]
    return sp.expand(R(n-1, k+1) - sum(
        (comb(n-2, j-1) * R(j, k) * R(n-j, 1)
         for j in range(1, n)), ZERO))


def exact_zero(expression):
    assert sp.Poly(sp.expand(expression), *xs, domain=sp.QQ).is_zero


def main():
    print('Python 3; SymPy ' + sp.__version__ + '; coefficient domain QQ')
    print('Source: OEIS A263634 revision 166, modified 2025-11-05; Kurkov 2025-03-30')
    L = logarithmic_polynomials()
    # EXAMPLE (1)--(6), transcribed from that revision, independently of R.
    source_examples = [
        'x1',
        '-x1**2+x2',
        '2*x1**3-3*x1*x2+x3',
        '-6*x1**4+12*x1**2*x2-4*x1*x3-3*x2**2+x4',
        '24*x1**5-60*x1**3*x2+20*x1**2*x3+30*x1*x2**2-5*x1*x4-10*x2*x3+x5',
        '-120*x1**6+360*x1**4*x2-120*x1**3*x3-270*x1**2*x2**2+30*x1**2*x4+120*x1*x2*x3+30*x2**3-6*x1*x5-15*x2*x4-10*x3**2+x6',
    ]
    source_rows = [
        [1], [-1, 1], [2, -3, 1], [-6, 12, -4, -3, 1],
        [24, -60, 20, 30, -5, -10, 1],
        [-120, 360, -120, -270, 30, 120, 30, -6, -15, -10, 1],
    ]
    local_names = {str(x): x for x in xs}
    for n, (text, row) in enumerate(zip(source_examples, source_rows), start=1):
        exact_zero(L[n] - sp.sympify(text, locals=local_names))
        coefficients = [int(c) for _, c in sp.Poly(L[n], *xs).terms(order='grlex')]
        assert coefficients == row, (n, coefficients, row)
        print('FIDELITY n=' + str(n) + ': EXAMPLE polynomial and coefficient row agree')
    for n in range(1, N+1):
        exact_zero(R(n, 1) - L[n])
        assert len(sp.Poly(L[n], *xs).terms()) == int(sp.partition(n))
        assert L[n].subs(dict.fromkeys(xs, ONE)) == (ONE if n == 1 else ZERO)
        print('SYMBOLIC_IDENTITY n=' + str(n) + ': R(n,1)-L_n = 0; terms=' +
              str(len(sp.Poly(L[n], *xs).terms())))
        print('L_' + str(n) + ' = ' + str(L[n]))
    # Independent finite jet of Q_k = (D^k F)/F, with F(0)=1.
    f = [ONE] + [xs[i-1]*sp.Rational(1, factorial(i)) for i in range(1, N+1)]
    inv = [ONE]
    for m in range(1, N+1):
        inv.append(sp.expand(-sum((f[i]*inv[m-i] for i in range(1, m+1)), ZERO)))
    count = 0
    for n in range(1, N+1):
        for k in range(1, N-n+2):
            # [t^(n-1)] D^k(F)/F uses source variables only through x_(n+k-1).
            q = sum((xs[k+a-1]*sp.Rational(1, factorial(a))*inv[n-1-a]
                     for a in range(n)), ZERO)
            exact_zero(factorial(n-1)*q - R(n, k))
            count += 1
    print('QUOTIENT_FAMILY: ' + str(count) + ' exact symbolic pairs; n,k >= 1, n+k-1 <= 8')
    print('Scope: finite exact checks only; the universal identity remains unproved here.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
