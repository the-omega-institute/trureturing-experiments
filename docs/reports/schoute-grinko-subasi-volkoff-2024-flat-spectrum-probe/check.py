#!/usr/bin/env python3
"""Checks for E. Schoute, D. Grinko, Y. Subasi and T. Volkoff, "Quantum programmable reflections", arXiv:2411.03648v2,
Section 4.2.1, Conjecture 4.6 / Eq. (183): for every n >= 1 there are probabilities q_j (j = (n mod 2)/2, ..., n/2) with
  (2J+1) / binom(n+2, 2) = sum_{j >= J/2} q_j / (2j+1) * | sum_{m=-j}^{j} C^{J0}_{jm, j -m} |^2
for every J in {(n mod 2), ..., n-2, n}.
Candidate explicit solution: q_{n/2-k} = (n-2k+1)^2 binom(2n-2k+1, n-k) binom(2k, k) / (4^n binom(n+2, 2)), k = 0..floor(n/2).
1. Exact (sympy Clebsch-Gordan, Condon-Shortley convention) for n = 1..12: the q are > 0, sum to 1, and satisfy every
   equation (183); the linear system has a unique solution.
2. Exact for 2j = 0..12: Lemma 1, U_{2j}(x) = sum_J |S_{j,J}|^2 P_J(x) with S_{j,J} = sum_m C^{J0}_{jm, j -m}.
3. Exact for n = 0..14: the positive connection formula G_n = C_n^{(3/2)} = sum_k a_{n,k} U_{n-2k}.
4. Exact: Racah's closed form in doubled integers (the convention preregistered in Lean, `clebschGordan a b c x y z`)
   equals sympy's Clebsch-Gordan coefficient for every 2j1, 2j2, 2J <= 6 and every admissible projection.
Prints ALL_OK.
"""
import sympy as sp
from sympy.physics.quantum.cg import CG
ok = True
x = sp.symbols('x')
half = sp.Rational(1, 2)
def S(j, J):
    return sp.nsimplify(sp.simplify(sum(CG(j, m, j, -m, J, 0).doit() for m in [j - i for i in range(int(2*j) + 1)])))
for n in range(1, 13):
    D = sp.binomial(n + 2, 2)
    js = [sp.Rational(n, 2) - k for k in range(n // 2 + 1)]
    q = {sp.Rational(n, 2) - k: (n - 2*k + 1)**2 * sp.binomial(2*n - 2*k + 1, n - k) * sp.binomial(2*k, k) / (4**n * D) for k in range(n // 2 + 1)}
    if any(v <= 0 for v in q.values()) or sum(q.values()) != 1: print('NORM_FAIL', n); ok = False
    Js = list(range(n % 2, n + 1, 2))
    rows = []
    for J in Js:
        coeffs = [ (S(j, J)**2 / (2*j + 1) if 2*j >= J else 0) for j in js ]
        rows.append(coeffs)
        lhs = sp.Rational(2*J + 1) / D; rhs = sum(c * q[j] for c, j in zip(coeffs, js))
        if sp.simplify(lhs - rhs) != 0: print('EQ_FAIL', n, J, lhs, rhs); ok = False
    M = sp.Matrix(rows)
    if M.rank() != len(js): print('NONUNIQUE', n, M.rank(), len(js)); ok = False
    print('n=%2d  D=%3d  q=%s' % (n, D, [str(q[j]) for j in js]), flush=True)
for tj in range(0, 13):
    j = sp.Rational(tj, 2)
    lhs = sp.chebyshevu(tj, x)
    rhs = sum(S(j, J)**2 * sp.legendre(J, x) for J in range(0, tj + 1))
    if sp.expand(lhs - rhs) != 0: print('LEMMA1_FAIL', j); ok = False
print('Lemma 1 checked for 2j = 0..12')
t = sp.symbols('t')
gen = sp.series((1 - 2*x*t + t**2)**sp.Rational(-3, 2), t, 0, 16).removeO()
for n in range(0, 15):
    Gn = sp.expand(gen.coeff(t, n))
    a = lambda k: (n - 2*k + 1) * sp.rf(sp.Rational(3, 2), n - k) * sp.rf(half, k) / (sp.factorial(k) * sp.factorial(n - k + 1))
    rhs = sum(a(k) * sp.chebyshevu(n - 2*k, x) for k in range(n // 2 + 1))
    if sp.expand(Gn - rhs) != 0 or any(a(k) <= 0 for k in range(n // 2 + 1)): print('LEMMA3_FAIL', n); ok = False
print('Lemma 3 (positive Gegenbauer-to-Chebyshev connection) checked for n = 0..14')

from math import factorial
def fz(v):
    return sp.Integer(factorial(max(v, 0)))
def racah(a, b, c, x, y, z):
    ok_ = (z == x + y and c <= a + b and a <= b + c and b <= a + c and (a + b + c) % 2 == 0 and abs(x) <= a and abs(y) <= b
           and abs(z) <= c and (a + x) % 2 == 0 and (b + y) % 2 == 0 and (c + z) % 2 == 0)
    if not ok_: return sp.Integer(0)
    pre = sp.sqrt((c + 1) * fz((c + a - b)//2) * fz((c - a + b)//2) * fz((a + b - c)//2) / fz((a + b + c)//2 + 1))
    pre *= sp.sqrt(fz((c + z)//2) * fz((c - z)//2) * fz((a - x)//2) * fz((a + x)//2) * fz((b - y)//2) * fz((b + y)//2))
    s = 0
    for k in range(a + b + c + 1):
        A = [(a + b - c)//2 - k, (a - x)//2 - k, (b + y)//2 - k, (c - b + x)//2 + k, (c - a - y)//2 + k]
        if min(A) >= 0:
            s += sp.Integer(-1)**k / (fz(k) * fz(A[0]) * fz(A[1]) * fz(A[2]) * fz(A[3]) * fz(A[4]))
    return sp.nsimplify(pre * s)
bad = 0; cnt = 0
for a in range(7):
    for b in range(7):
        for c in range(7):
            for x in range(-a, a + 1, 2):
                for y in range(-b, b + 1, 2):
                    z = x + y
                    if abs(z) > c: continue
                    cnt += 1
                    ref = sp.nsimplify(sp.simplify(CG(sp.Rational(a, 2), sp.Rational(x, 2), sp.Rational(b, 2), sp.Rational(y, 2), sp.Rational(c, 2), sp.Rational(z, 2)).doit()))
                    if sp.simplify(racah(a, b, c, x, y, z) - ref) != 0: bad += 1
print('Racah doubled-integer convention vs sympy CG: %d cases, %d mismatches' % (cnt, bad))
if bad: ok = False
print('ALL_OK' if ok else 'FAIL')
