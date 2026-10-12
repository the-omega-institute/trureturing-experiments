#!/usr/bin/env python3
"""OEIS A308815 (Mathar's conjecture, Jan 27 2020): with C(m) = binomial(2m,m)/(m+1) and
a(n) = 16^n * sum_{i=0..n} (-2)^i binomial(n,i) C(2n-i)   (the entry's data and PARI program; %N prints 2^i),
n(2n+1)a(n) - 32(9n-4)(n-1)a(n-1) + 1024(3n^2-23n+29)a(n-2) + 65536(n-2)(2n-5)a(n-3) = 0 for n >= 3.
Exact checks (not a proof):
  A. the (-2)^i definition reproduces the entry's 17 data terms; the printed 2^i definition does not;
  B. the recurrence holds exactly for 3 <= n <= 300;
  C. B(z) = sum b_n z^n, b_n = a_n/16^n, satisfies 4z(1+z)B^4 - (1+4z)B^2 + 1 = 0 to order z^120;
  D. the differential operator L (theta = z d/dz) satisfies L(B) = 2Bz(2z-1)^3 (8B^2z^2+8B^2z-12z-3) P(z,B)/D^3
     as a rational-function identity once theta B is replaced by its implicit-differentiation value.
Exit 0 and final line `bad= 0` when every check agrees."""
import sys
from math import comb
import sympy as sp

bad = 0
def cat(m): return comb(2 * m, m) // (m + 1)
def a(n, sgn=-2): return 16 ** n * sum(sgn ** i * comb(n, i) * cat(2 * n - i) for i in range(n + 1))
data = [1, 0, 512, 32768, 3014656, 285212672, 28521267200, 2950642532352, 313455303196672, 33990302461067264,
        3747096045042008064, 418698249981751459840, 47315379509207945445376, 5398203160147467800936448,
        620943817223351793348509696, 71934849808689842265157271552, 8385447350799836001482047488000]
okA = [a(n) for n in range(17)] == data and [a(n, 2) for n in range(17)] != data
print('A: (-2)^i matches data:', [a(n) for n in range(17)] == data, '; printed 2^i matches:', [a(n, 2) for n in range(17)] == data)
if not okA: bad += 1

A = [a(n) for n in range(301)]
fails = [n for n in range(3, 301) if n * (2 * n + 1) * A[n] - 32 * (9 * n - 4) * (n - 1) * A[n - 1]
         + 1024 * (3 * n * n - 23 * n + 29) * A[n - 2] + 65536 * (n - 2) * (2 * n - 5) * A[n - 3] != 0]
print('B: recurrence failures for 3<=n<=300:', fails[:5], 'count', len(fails))
if fails: bad += 1

z = sp.symbols('z')
N = 120
b = [sp.Rational(A[n], 16 ** n) if n <= 300 else None for n in range(N + 1)]
Bs = sum(b[n] * z ** n for n in range(N + 1))
P = sp.expand(4 * z * (1 + z) * Bs ** 4 - (1 + 4 * z) * Bs ** 2 + 1)
poly = sp.Poly(P, z)
low = [poly.coeff_monomial(z ** k) for k in range(N + 1)]
print('C: nonzero coefficients of P(z,B) below z^%d:' % (N + 1), sum(1 for c in low if c != 0))
if any(c != 0 for c in low): bad += 1

Bv, zv = sp.symbols('B z')
Pz = 4 * zv * (1 + zv) * Bv ** 4 - (1 + 4 * zv) * Bv ** 2 + 1
Dz = 8 * Bv ** 2 * zv ** 2 + 8 * Bv ** 2 * zv - 4 * zv - 1
thetaB = -2 * Bv * zv * (2 * Bv ** 2 * zv + Bv ** 2 - 1) / Dz          # implicit differentiation value
# check thetaB against -z P_z / P_B on the curve P = 0 (both are theta B there): difference is a multiple of P
impl = -zv * sp.diff(Pz, zv) / sp.diff(Pz, Bv)
r1 = sp.simplify(sp.together(impl - thetaB))
num1 = sp.numer(sp.together(r1))
q1, rem1 = sp.div(sp.Poly(num1, Bv, zv), sp.Poly(Pz, Bv, zv))
# theta applied to a function f(B,z): theta f = z f_z + f_B * thetaB
def theta(f): return zv * sp.diff(f, zv) + sp.diff(f, Bv) * thetaB
t1 = thetaB
t2 = sp.together(theta(t1))
Lexpr = (2 * t2 + t1) + zv * (-18 * t2 - 10 * t1) + zv ** 2 * (12 * t2 - 44 * t1 - 20 * Bv) + zv ** 3 * (32 * t2 + 48 * t1 + 16 * Bv)
claimed = 2 * Bv * zv * (2 * zv - 1) ** 3 * (8 * Bv ** 2 * zv ** 2 + 8 * Bv ** 2 * zv - 12 * zv - 3) * Pz / Dz ** 3
diff = sp.simplify(sp.together(Lexpr - claimed))
print('D: thetaB equals -zP_z/P_B modulo P:', rem1.is_zero, '; L(B) minus the claimed factorized form:', diff)
if not (rem1.is_zero and diff == 0): bad += 1
print('bad=', bad)
sys.exit(0 if bad == 0 else 1)
