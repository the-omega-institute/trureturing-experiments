"""Check of the proposed refutation of Zhang-Xiang-Li-Jost-Fei (arXiv:1912.12321, PRE 100 (2019) 062139),
Section III.C: 'The incompatibility probability for a pair of general qubit measurements is conjectured to be 1/4'.

(1) Exact reduction: P = 9 I with I = 4 int_[0,1]^4 rho^3 sqrt(1-rho) a(1-a)(1-t^2)(1-s^2) sqrt(1-rho(a t^2+(1-a) s^2)).
    I is computed by tensor Gauss-Legendre quadrature (rho through z = sqrt(1-rho) to remove the endpoint singularity).
(2) Rigorous bound: P <= 5910272/23648625 = 1/4 - 7537/94594500 (exact rational arithmetic of the moment certificate).
(3) Monte Carlo of the source's definition: two binary qubit POVMs E(x,a) sampled uniformly from |x|+|a|<=1,
    incompatible by the Yu-Liu-Li-Oh criterion (arXiv:0805.1538 Thm 1); the estimate must agree with 9 I.
"""
import numpy as np
from fractions import Fraction as Fr
bad = 0
# (2) exact certificate arithmetic
B0, B1, B2 = Fr(32, 315), Fr(256, 3465), Fr(512, 9009)
C0, C1, C2 = Fr(2, 27), Fr(2, 135), Fr(118, 23625)
# moments recomputed exactly
def mom_t(k): return Fr(2, (2*k+1)*(2*k+3))   # int_0^1 t^{2k}(1-t^2)
def mom_a(i, j):                               # int_0^1 a^i (1-a)^j
    from math import factorial
    return Fr(factorial(i)*factorial(j), factorial(i+j+1))
c0 = mom_a(1,1)*mom_t(0)**2
c1 = mom_a(2,1)*mom_t(1)*mom_t(0) + mom_a(1,2)*mom_t(0)*mom_t(1)
c2 = mom_a(3,1)*mom_t(2)*mom_t(0) + 2*mom_a(2,2)*mom_t(1)*mom_t(1) + mom_a(1,3)*mom_t(0)*mom_t(2)
def bj(j):   # int_0^1 rho^{3+j} sqrt(1-rho): z = sqrt(1-rho), rho = 1-z^2, drho = -2z dz -> int_0^1 (1-z^2)^{3+j} z * 2z dz
    from math import comb
    n = 3 + j
    return sum(Fr(comb(n, k) * (-1)**k * 2, 2*k + 3) for k in range(n + 1))
if (c0, c1, c2) != (C0, C1, C2): bad += 1; print("FAIL C moments", c0, c1, c2)
if (bj(0), bj(1), bj(2)) != (B0, B1, B2): bad += 1; print("FAIL B moments", bj(0), bj(1), bj(2))
Ibound = 4 * (B0*C0 - B1*C1/2 - B2*C2/8)
Pbound = 9 * Ibound
print("bound P <=", Pbound, "=", float(Pbound), "; 1/4 - bound =", Fr(1,4) - Pbound)
if Pbound != Fr(5910272, 23648625) or Fr(1,4) - Pbound != Fr(7537, 94594500): bad += 1; print("FAIL bound arithmetic")
# sqrt(1-z) <= p(z) on [0,1]
z = np.linspace(0, 1, 100001)
if np.any(np.sqrt(1 - z) > 1 - z/2 - z**2/8 + 1e-15): bad += 1; print("FAIL p(z) certificate")
# (1) quadrature of I
def gl(n):
    x, w = np.polynomial.legendre.leggauss(n); return (x + 1) / 2, w / 2
n = 48
zz, wz = gl(n); aa, wa = gl(n); tt, wt = gl(n); ss, ws = gl(n)
rho = 1 - zz**2; wrho = wz * 2 * zz          # drho = 2 z dz ; sqrt(1-rho) = z
f_rho = rho**3 * zz * wrho
A = aa[:, None, None]; T = tt[None, :, None]; Sx = ss[None, None, :]
W = A * T**2 + (1 - A) * Sx**2
wfac = (A * (1 - A) * (1 - T**2) * (1 - Sx**2)) * (wa[:, None, None] * wt[None, :, None] * ws[None, None, :])
I = 0.0
for r, fr in zip(rho, f_rho):
    I += fr * np.sum(wfac * np.sqrt(np.clip(1 - r * W, 0, None)))
I *= 4
P = 9 * I
print("quadrature P = 9 I =", P, "; 1/4 - P =", 0.25 - P)
if not (P <= float(Pbound) + 1e-12 and P < 0.25): bad += 1; print("FAIL quadrature above bound")
# (3) Monte Carlo of the source definition
rng = np.random.default_rng(1912)
def sample(N):
    # uniform on {(x,a): |x|+|a|<=1} in R x R^3 by rejection from the box
    out = np.empty((0, 4))
    while len(out) < N:
        c = rng.uniform(-1, 1, size=(2 * N, 4))
        keep = np.abs(c[:, 0]) + np.linalg.norm(c[:, 1:], axis=1) <= 1
        out = np.vstack([out, c[keep]])
    return out[:N]
def F(x, r):
    return (np.sqrt(np.clip((1 + x)**2 - r**2, 0, None)) + np.sqrt(np.clip((1 - x)**2 - r**2, 0, None))) / 2
N = 4_000_000; inc = 0; tot = 0
for _ in range(5):
    m1 = sample(N); m2 = sample(N)
    x, a = m1[:, 0], m1[:, 1:]; y, b = m2[:, 0], m2[:, 1:]
    Fx = F(x, np.linalg.norm(a, axis=1)); Fy = F(y, np.linalg.norm(b, axis=1))
    lhs = (1 - Fx**2 - Fy**2) * (1 - x**2 / Fx**2 - y**2 / Fy**2)
    rhs = (np.sum(a * b, axis=1) - x * y)**2
    inc += np.count_nonzero(lhs > rhs); tot += N
Pmc = inc / tot; se = np.sqrt(Pmc * (1 - Pmc) / tot)
print(f"Monte Carlo P = {Pmc:.5f} +- {se:.5f} over {tot} pairs; quadrature {P:.5f}")
if abs(Pmc - P) > 4 * se: bad += 1; print("FAIL Monte Carlo disagrees with the reduction")
print("bad=", bad)
raise SystemExit(1 if bad else 0)
