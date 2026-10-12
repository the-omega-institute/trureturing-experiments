#!/usr/bin/env python3
"""Wang-Chen-Song, Phys. Scr. 101 (2026) 315210 (arXiv:2603.20717v1), Conclusions, fourth question: the achievable upper
bound of the number of distinct eigenvalues of extreme two-qutrit absolutely PPT states. Proposed answer 6.
Checks (not a proof):
  A. exact, for the family a..f(h), h in [0, 1/100]: six strictly decreasing levels; L1(lambda) p = 0 and
     L2(lambda) q(h) = 0 with p = (13,11,8), q(h) = (52, 39+13h, 34+22h); both pencils PSD of rank two;
     det M(h) = -512 Z(h) != 0;
  B. exact: w . mu = p^T L1(mu) p for every spectrum mu, with w = (-286,-208,-176,128,176,208,242,286,338) increasing;
  C. numerical sanity of the spectral criterion: at h = 0 no sampled or locally optimised global unitary makes the
     partial transpose of U rho U^* negative; a spectrum violating L1 >= 0 is caught by the same search.
Exit 0 and final line `bad= 0` when every check agrees."""
import sys
import numpy as np
import sympy as sp

bad = 0
h = sp.symbols('h')
a = 3498 + 20188*h + 7706*h**2; b = 3432 + 4264*h - 624*h**2; c = 2200 - 2200*h - 3168*h**2
d = 1859 - 6760*h - 507*h**2; e = 1628 + 5808*h - 1452*h**2; f = 1430 + 5668*h + 2470*h**2
Z = sp.expand(a + b + c + d + 2*e + 3*f)
lam = [a, b, c, d, e, e, f, f, f]                         # unnormalized, decreasing
def L1(l): return sp.Matrix([[2*l[8], l[7]-l[0], l[5]-l[1]], [l[7]-l[0], 2*l[6], l[4]-l[2]], [l[5]-l[1], l[4]-l[2], 2*l[3]]])
def L2(l): return sp.Matrix([[2*l[8], l[7]-l[0], l[6]-l[1]], [l[7]-l[0], 2*l[5], l[4]-l[2]], [l[6]-l[1], l[4]-l[2], 2*l[3]]])
p = sp.Matrix([13, 11, 8]); q = sp.Matrix([52, 39 + 13*h, 34 + 22*h])
A_, B_ = L1(lam), L2(lam)
okA = sp.expand(Z - (18535 + 44112*h + 7913*h**2)) == 0
okA &= all(sp.expand(x) == 0 for x in A_ * p) and all(sp.expand(x) == 0 for x in B_ * q)
okA &= sp.expand(A_.det()) == 0 and sp.expand(B_.det()) == 0
gaps = [sp.expand(a-b), sp.expand(b-c), sp.expand(c-d), sp.expand(d-e), sp.expand(e-f), f]
minorA = sp.expand(A_[:2, :2].det()); minorB = sp.expand(B_[:2, :2].det())
def positive_on_interval(poly):
    P = sp.Poly(poly, h)
    roots = [r for r in P.real_roots() if 0 <= r <= sp.Rational(1, 100)]
    return not roots and poly.subs(h, 0) > 0 and poly.subs(h, sp.Rational(1, 100)) > 0
for g in gaps + [minorA, minorB, A_[0, 0], B_[0, 0]]:
    okA &= positive_on_interval(g)
x = sp.symbols('x1:7')
t = [x[0], x[1], x[2], x[3], x[4], x[4], x[5], x[5], x[5]]
rows = list(L1(t) * p) + list((L2(t) * q)[:2]) + [sum(t)]
M = sp.Matrix([[sp.expand(r).coeff(xi) for xi in x] for r in rows])
detM = sp.expand(M.det())
okA &= sp.expand(detM - (-512) * Z) == 0 and positive_on_interval(Z)
print('A: family identities, PSD rank-two pencils, gaps, det M(h) = -512 Z(h):', 'ok' if okA else 'BAD')
if not okA: bad += 1

mu = sp.symbols('m1:10')
w = [-286, -208, -176, 128, 176, 208, 242, 286, 338]
okB = sp.expand(sum(wi*mi for wi, mi in zip(w, mu)) - (p.T * L1(list(mu)) * p)[0]) == 0 and all(w[i] < w[i+1] for i in range(8))
print('B: w . mu = p^T L1(mu) p and w increasing:', 'ok' if okB else 'BAD')
if not okB: bad += 1

rng = np.random.default_rng(20261012)
def pt_min_eig(rho, U):
    R = U @ rho @ U.conj().T
    R4 = R.reshape(3, 3, 3, 3)
    PT = R4.transpose(0, 3, 2, 1).reshape(9, 9)            # transpose on the second factor
    return np.linalg.eigvalsh((PT + PT.conj().T) / 2).min()
def haar(n):
    Zm = (rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))) / np.sqrt(2)
    Q, R = np.linalg.qr(Zm); return Q * (np.diag(R) / np.abs(np.diag(R)))
def search(spec, trials=3000, steps=60):
    rho = np.diag(np.array(spec, dtype=float))
    best = min(pt_min_eig(rho, haar(9)) for _ in range(trials))
    for _ in range(40):                                      # crude local descent from random starts
        U = haar(9); val = pt_min_eig(rho, U); step = 0.3
        for _ in range(steps):
            Hm = rng.normal(size=(9, 9)) + 1j * rng.normal(size=(9, 9)); Hm = (Hm + Hm.conj().T) / 2
            w_, V = np.linalg.eigh(Hm); Uc = V @ np.diag(np.exp(1j * step * w_)) @ V.conj().T @ U
            vc = pt_min_eig(rho, Uc)
            if vc < val: U, val = Uc, vc
            else: step *= 0.9
        best = min(best, val)
    return best
spec0 = [float(v) / 1685 for v in [318, 312, 200, 169, 148, 148, 130, 130, 130]]
m0 = search(spec0)
bad_spec = [0.30, 0.25, 0.15, 0.1, 0.05, 0.05, 0.04, 0.03, 0.03]; s = sum(bad_spec); bad_spec = [v / s for v in bad_spec]
L1bad = np.array([[2*bad_spec[8], bad_spec[7]-bad_spec[0], bad_spec[5]-bad_spec[1]], [bad_spec[7]-bad_spec[0], 2*bad_spec[6], bad_spec[4]-bad_spec[2]], [bad_spec[5]-bad_spec[1], bad_spec[4]-bad_spec[2], 2*bad_spec[3]]])
m1 = search(bad_spec, trials=1500, steps=40)
okC = m0 > -1e-10 and np.linalg.eigvalsh(L1bad).min() < 0 and m1 < -1e-4
print(f'C: min PT eigenvalue found at h=0: {m0:.3e} (expect >= 0); control spectrum with L1 min eig {np.linalg.eigvalsh(L1bad).min():.3e}: found {m1:.3e} (expect < 0)', 'ok' if okC else 'BAD')
if not okC: bad += 1
print('bad=', bad)
sys.exit(0 if bad == 0 else 1)
