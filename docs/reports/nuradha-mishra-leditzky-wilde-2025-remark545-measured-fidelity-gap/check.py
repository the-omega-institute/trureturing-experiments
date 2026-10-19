"""Nuradha-Mishra-Leditzky-Wilde, arXiv:2404.16101v3, Remark 5.45: is there a tuple with F_M > F_U?

Witness: R = I/2, A = (I + Z/2)/2, B = (I + X/2)/2 (qubits).
1. F_U from the trace-norm definition F(rho, sigma) = ||sqrt(rho) sqrt(sigma)||_1 (numerically) and in closed
   form (2 sqrt 6 + 2 sqrt 2 + sqrt 14)/12.
2. F_M upper estimate: LP over POVMs with rank-one effects on N random Bloch directions (inf over a subset).
3. F_M lower bound by a rational dual certificate: for every effect E >= 0,
   g(E) = sum_{i<j} sqrt(tr(E rho_i) tr(E rho_j)) >= tr(E Y) with Y = y0 I + y (X + Z)/... ;
   g is concave and 1-homogeneous, so it suffices to check rank-one effects E = (I + n.sigma)/2, i.e.
   h(n) >= y0 + y1 (n_x + n_z) on the unit sphere; then sum_x g(E_x) >= tr(Y) = 2 y0 for every POVM.
   The check uses a fine (theta, phi) grid plus a Lipschitz margin.
Requires numpy, scipy.  Usage: python3 check.py  (prints RESULT PASS / RESULT FAIL).
"""
import numpy as np
from fractions import Fraction as Fr
from scipy.linalg import sqrtm
from scipy.optimize import linprog

I = np.eye(2); X = np.array([[0, 1], [1, 0]]); Y = np.array([[0, -1j], [1j, 0]]); Z = np.diag([1, -1])
R = I / 2; A = (I + Z / 2) / 2; B = (I + X / 2) / 2
def F(r, s):
    return float(np.sum(np.linalg.svd(sqrtm(r) @ sqrtm(s), compute_uv=False)).real)

def h(n):
    pr = 0.5; pa = (1 + n[..., 2] / 2) / 2; pb = (1 + n[..., 0] / 2) / 2
    return np.sqrt(pr * pa) + np.sqrt(pr * pb) + np.sqrt(pa * pb)

def main():
    ok = True
    FU = (F(R, A) + F(R, B) + F(A, B)) / 3
    FU_closed = (2 * np.sqrt(6) + 2 * np.sqrt(2) + np.sqrt(14)) / 12
    print(f"F_U = {FU:.12f} (trace norm), closed form {FU_closed:.12f}")
    ok &= abs(FU - FU_closed) < 1e-10
    rng = np.random.default_rng(1)
    v = rng.normal(size=(20000, 3)); v /= np.linalg.norm(v, axis=1)[:, None]
    v = np.vstack([v, -v, np.eye(3), -np.eye(3)])
    res = linprog(h(v), A_eq=np.vstack([np.ones(len(v)), v.T]), b_eq=[2, 0, 0, 0], bounds=(0, None), method="highs")
    FM_up = res.fun / 3
    print(f"F_M <= {FM_up:.12f} (LP over 40006 rank-one directions); LP dual y0 = {res.eqlin.marginals[0]:.6f}")
    y0, y1 = Fr(36, 25), Fr(53, 200)
    th = np.linspace(0, np.pi, 2001); ph = np.linspace(0, 2 * np.pi, 4001)
    T, P = np.meshgrid(th, ph)
    n = np.stack([np.sin(T) * np.cos(P), np.sin(T) * np.sin(P), np.cos(T)], -1)
    slack = h(n) - (float(y0) + float(y1) * (n[..., 0] + n[..., 2]))
    # every unit vector lies within geodesic distance sqrt(dtheta^2 + dphi^2)/2 of a grid point, and
    # h - linear is Lipschitz on the sphere with constant <= 3/4 + 2 y1 (each square-root term has gradient <= 1/4)
    lip = 0.75 + 2 * float(y1); margin = lip * np.hypot(np.pi / 2000, 2 * np.pi / 4000) / 2
    print(f"certificate y0 = {y0}, y1 = {y1}: min slack on grid {slack.min():.6e}, Lipschitz margin {margin:.6e}")
    ok &= slack.min() > margin
    FM_low = float(2 * y0) / 3
    print(f"F_M >= 2 y0 / 3 = {FM_low:.12f} > F_U = {FU:.12f}: {FM_low > FU}; gap >= {FM_low - FU:.6e}")
    ok &= FM_low > FU and FM_low <= FM_up + 1e-9
    print("RESULT", "PASS" if ok else "FAIL")

if __name__ == "__main__":
    main()
