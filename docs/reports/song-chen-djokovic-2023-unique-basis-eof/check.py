#!/usr/bin/env python3
"""Exact checks for Song-Chen-Djokovic, arXiv:2307.01980v3 (IEEE TIT 70 (2024) 8806), p.2-3: G (Eq. (2)), Psi_1..3
(Eq. (3)) and the question 'We don't know whether this EOF is exactly equal to that of the pure state |Psi_1>.'

1. Psi_j (built from the columns r_i of G) are orthonormal; G^dagger G is a dd-matrix with (G^dagger G)_B = I_3.
2. Tr_B |Psi_j><Psi_j| = cyclic shifts of diag(4/9, 1/3, 2/9); H0 = H(4/9, 1/3, 2/9) bits.
3. phi_+- = (Psi_1 +- sqrt(3) Psi_2)/2: unit norm, <phi_+|phi_-> = -1/2, (phi_+ phi_+^dag + phi_- phi_-^dag)/2 = (P1 + 3 P2)/4,
   marginal spectra exactly (17, 10, 9)/36.
4. H* = H(17/36, 10/36, 9/36) < H0 by the integer inequality 17^17 5^10 3^6 > 2^102.
5. rho0 = (P1 + P2 + P3)/3 = (4/9)(P1 + 3P2)/4 + (P1 (2/9) + P2 (0) + P3 (1/3)) ... explicit ensemble with cost
   (4/9) H* + (5/9) H0 < H0.
Prints ALL_OK.
"""
import sympy as sp
I = sp.I; s3 = sp.sqrt(3); s2 = sp.sqrt(2); s11 = sp.sqrt(11); s33 = sp.sqrt(33)
G = sp.zeros(9, 9)
G[0,0] = sp.Rational(2,3); G[1,1] = 1/s3; G[1,5] = I/(3*s3); G[1,6] = 2/(3*s3); G[2,2] = s2/3; G[3,3] = 1/s3
G[4,4] = s2/3; G[5,5] = s11/(3*s3); G[5,6] = (3+2*I)/(3*s33); G[6,6] = 1/s33; G[7,7] = sp.Rational(2,3); G[8,8] = 1/s3
r = [G[:, i] for i in range(9)]
def psi(cols):  # |1>⊗r_a + |2>⊗r_b + |3>⊗r_c  -> vector indexed (a, b) -> 9a + b
    v = sp.zeros(27, 1)
    for a, c in enumerate(cols):
        for b in range(9): v[9*a + b] = r[c][b]
    return v
Psi = [psi([0, 3, 6]), psi([1, 4, 7]), psi([2, 5, 8])]
ok = True
def ip(u, v): return sp.simplify((u.H * v)[0])
for i in range(3):
    for j in range(3):
        if sp.simplify(ip(Psi[i], Psi[j]) - (1 if i == j else 0)) != 0: print('NOT_ORTHONORMAL', i, j); ok = False
GG = (G.H * G).applyfunc(sp.simplify)
MB = sum((GG[3*k:3*k+3, 3*k:3*k+3] for k in range(3)), sp.zeros(3, 3))
if MB != sp.eye(3): print('MB_NOT_I3'); ok = False
def trB(v):
    R = sp.zeros(3, 3)
    for a in range(3):
        for a2 in range(3):
            R[a, a2] = sp.simplify(sum(v[9*a + b] * sp.conjugate(v[9*a2 + b]) for b in range(9)))
    return R
targets = [sp.diag(sp.Rational(4,9), sp.Rational(1,3), sp.Rational(2,9)), sp.diag(sp.Rational(1,3), sp.Rational(2,9), sp.Rational(4,9)),
           sp.diag(sp.Rational(2,9), sp.Rational(4,9), sp.Rational(1,3))]
for j in range(3):
    if (trB(Psi[j]) - targets[j]).applyfunc(sp.simplify) != sp.zeros(3, 3): print('MARGINAL', j, trB(Psi[j])); ok = False
php = (Psi[0] + s3 * Psi[1]) / 2; phm = (Psi[0] - s3 * Psi[1]) / 2
if sp.simplify(ip(php, php) - 1) != 0 or sp.simplify(ip(phm, phm) - 1) != 0: print('PHI_NORM'); ok = False
if sp.simplify(ip(php, phm) + sp.Rational(1, 2)) != 0: print('PHI_OVERLAP'); ok = False
P = [v * v.H for v in Psi]
mix = ((php * php.H + phm * phm.H) / 2 - (P[0] + 3 * P[1]) / 4).applyfunc(sp.simplify)
if mix != sp.zeros(27, 27): print('ENSEMBLE_MIX'); ok = False
lam = sp.symbols('lam')
for ph, nm in ((php, '+'), (phm, '-')):
    R = trB(ph)
    cp = sp.factor(sp.expand((36 * R - lam * sp.eye(3)).det()))
    target = sp.expand((17 - lam) * (10 - lam) * (9 - lam))
    if sp.expand(cp) != target: print('SPECTRUM', nm, cp); ok = False
L = 17**17 * 5**10 * 3**6; Rr = 2**102
if not L > Rr: print('ENTROPY_INEQ'); ok = False
print('L - R =', L - Rr)
H = lambda ps: -sum(p * sp.log(p, 2) for p in ps)
H0 = H([sp.Rational(4,9), sp.Rational(1,3), sp.Rational(2,9)]); Hs = H([sp.Rational(17,36), sp.Rational(10,36), sp.Rational(9,36)])
print('H0 = %.10f  H* = %.10f' % (float(H0), float(Hs)))
# rho0 = (P1+P2+P3)/3 = (4/9) * (P1 + 3 P2)/4 + (2/9) P1 + (1/3) P3 ; P2 weight: (4/9)(3/4) = 1/3 exactly
rho0 = (P[0] + P[1] + P[2]) / 3
dec = (sp.Rational(4,9) * (P[0] + 3*P[1]) / 4 + sp.Rational(2,9) * P[0] + sp.Rational(1,3) * P[2] - rho0).applyfunc(sp.simplify)
if dec != sp.zeros(27, 27): print('RHO0_DECOMP'); ok = False
cost = sp.Rational(4,9) * Hs + sp.Rational(5,9) * H0
print('cost(rho0 ensemble) = %.10f < H0 = %.10f' % (float(cost), float(H0)))
if not float(cost) < float(H0): ok = False
print('ALL_OK' if ok else 'FAIL')
