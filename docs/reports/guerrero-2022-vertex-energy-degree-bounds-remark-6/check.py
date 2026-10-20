"""Guerrero, arXiv:2201.01252v1, Remark 6 (Conjecture): d_min LE_G(v) <= E_G(v) <= d_max LE_G(v) fails.

E_G(v) = |A|_vv and LE_G(v) = |L - I|_vv with L = I - D^{-1/2} A D^{-1/2}. Checks:
1. numerically, for the path P6 and its blow-ups G_t (each vertex replaced by t independent copies, each edge by K_{t,t}),
   t = 1..4, that the central vertices violate E_G(v) <= d_max LE_G(v);
2. exactly (sympy), the closed forms E_{P6}(3) = (7 + 4b - 2b^2)/7 with b = 2cos(2pi/7) and LE_{P6}(3) = (4 + sqrt 5)/10,
   and the chain E_{P6}(3) > 79/63 > 5/4 > 2 LE_{P6}(3).
"""
import sys
import numpy as np
import sympy as sp


def absm(M):
    w, V = np.linalg.eigh(M)
    return V @ np.diag(np.abs(w)) @ V.T


def energies(A):
    d = A.sum(1)
    Dm = np.diag(d ** -0.5)
    R = Dm @ A @ Dm
    return np.diag(absm(A)), np.diag(absm(R)), d.min(), d.max()


P = np.zeros((6, 6))
for i in range(5):
    P[i, i + 1] = P[i + 1, i] = 1
ok = True
for t in (1, 2, 3, 4):
    A = np.kron(P, np.ones((t, t)))
    E, LE, dmin, dmax = energies(A)
    bad = [v for v in range(6 * t) if E[v] > dmax * LE[v] + 1e-12]
    centre = [2 * t + j for j in range(2 * t)]
    hit = sorted(bad) == centre
    ok = ok and hit
    print(f"t={t}: n={6 * t}, d_min={dmin:.0f}, d_max={dmax:.0f}, upper-bound violations at the 2t central vertices only: {hit}; "
          f"E={E[2 * t]:.6f}, d_max*LE={dmax * LE[2 * t]:.6f}")
b = 2 * sp.cos(2 * sp.pi / 7)
e3_closed = (7 + 4 * b - 2 * b ** 2) / 7
r3_closed = (4 + sp.sqrt(5)) / 10
E_num, LE_num, _, _ = energies(P)
close_e = abs(float(e3_closed) - E_num[2]) < 1e-12
close_r = abs(float(r3_closed) - LE_num[2]) < 1e-12
chain = float(e3_closed) > 79 / 63 > 5 / 4 > float(2 * r3_closed)
ok = ok and close_e and close_r and chain
print(f"closed forms match numerics: E(3) {close_e}, LE(3) {close_r}; E(3) = {float(e3_closed):.12f} > 79/63 = {79 / 63:.12f} > 5/4 > 2 LE(3) = {float(2 * r3_closed):.12f}: {chain}")
print(f"b = 2cos(2pi/7) satisfies b^3 + b^2 - 2b - 1 = 0 numerically: {abs(float(b ** 3 + b ** 2 - 2 * b - 1)) < 1e-12}")
print("RESULT PASS" if ok else "RESULT FAIL")
sys.exit(0 if ok else 1)
