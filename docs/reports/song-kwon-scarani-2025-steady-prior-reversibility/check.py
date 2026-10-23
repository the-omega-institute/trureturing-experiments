#!/usr/bin/env python3
"""Checks for Song, Kwon and Scarani, arXiv:2510.26895v1, Section V open question 2 (numpy only).

N(X) = Tr_E[U (X (x) xi) U^dag] (1); Petz map (2) with N^dag(Y) = Tr_E[U^dag (Y (x) 1) U (1 (x) xi)]; tabletop reverse
map Tr_E[U^dag (X (x) xi') U] (3). TTR (4) needs some xi' with Petz = tabletop reverse; here the least-squares residual
over all complex 3x3 (2x2) matrices xi' is computed, so a positive residual excludes every xi', density or not.
Fidelity: Examples 1 and 2 of the source are TTR (residual ~ 0). Claim: the family S = E = C^3, xi = diag(1/2,1/3,1/6),
gamma = diag(17,13,14)/44, U_t = exp(-itH) keeps gamma steady for all t and is not TTR for t in (0, pi/2]; at
t = pi/2 the Petz population matrix equals the exact rationals of the prereg.
Exits 1 on failure.
"""
import sys
import numpy as np

def ptrE(M, ds, de):
    return np.einsum('aebe->ab', M.reshape(ds, de, ds, de))

def hsqrt(P):
    w, V = np.linalg.eigh((P + P.conj().T) / 2)
    return V @ np.diag(np.sqrt(np.clip(w, 0, None))) @ V.conj().T

def residual(U, xi, g):
    ds, de = g.shape[0], xi.shape[0]
    N = lambda X: ptrE(U @ np.kron(X, xi) @ U.conj().T, ds, de)
    Nadj = lambda Y: ptrE(U.conj().T @ np.kron(Y, np.eye(de)) @ U @ np.kron(np.eye(ds), xi), ds, de)
    Ng = N(g); r = np.linalg.inv(hsqrt(Ng)); sg = hsqrt(g)
    basis = lambda n: [np.eye(n)[:, [i]] @ np.eye(n)[[j], :] for i in range(n) for j in range(n)]
    A, b = [], []
    for X in basis(ds):
        b.append((sg @ Nadj(r @ X @ r) @ sg).ravel())
        A.append(np.array([ptrE(U.conj().T @ np.kron(X, Ep) @ U, ds, de).ravel() for Ep in basis(de)]).T)
    A = np.vstack(A); b = np.concatenate(b)
    sol, *_ = np.linalg.lstsq(A, b, rcond=None)
    return float(np.linalg.norm(A @ sol - b)), float(np.linalg.norm(Ng - g))

def unitary_exp(H, t):
    w, V = np.linalg.eigh(H)
    return V @ np.diag(np.exp(-1j * t * w)) @ V.conj().T

ok = True
rng = np.random.default_rng(20261023)
def rand_herm(n):
    A = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)); return (A + A.conj().T) / 2
P = [np.diag(v) for v in ([1., 1., 0.], [0., 0., 1.])]
A = np.zeros((3, 3), complex); A[:2, :2] = rand_herm(2)
H1 = np.kron(P[0], rand_herm(3)) + np.kron(P[1], rand_herm(3)) + np.kron(A, rand_herm(3))
g1 = np.diag([0.3, 0.3, 0.4])
r1, s1 = residual(unitary_exp(H1, 1.0), np.eye(3) / 3, g1)
print(f"fidelity Example 1 (unital, xi = 1/3, steady gamma = diag(0.3, 0.3, 0.4), U commuting with gamma (x) 1): steady {s1:.1e}, TTR residual {r1:.1e}")
ok &= r1 < 1e-10 and s1 < 1e-10
CNOT = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [0, 0, 1, 0]], complex)
r2, s2 = residual(CNOT, np.diag([0.7, 0.3]), np.diag([0.6, 0.4]))
print(f"fidelity Example 2 (controlled-X, diagonal xi and gamma): steady {s2:.1e}, TTR residual {r2:.1e}")
ok &= r2 < 1e-10 and s2 < 1e-10
H = np.zeros((9, 9))
for a, b in (((0, 0), (1, 1)), ((1, 0), (2, 2)), ((2, 0), (0, 2))):
    i, j = 3 * a[0] + a[1], 3 * b[0] + b[1]; H[i, j] = H[j, i] = 1
xi = np.diag([1 / 2, 1 / 3, 1 / 6]); g = np.diag([17, 13, 14]) / 44
worst_steady, min_res = 0.0, np.inf
for t in np.linspace(np.pi / 100, np.pi / 2, 50):
    r, s = residual(unitary_exp(H, t), xi, g)
    worst_steady = max(worst_steady, s); min_res = min(min_res, r)
r_half, _ = residual(unitary_exp(H, np.pi / 2), xi, g)
print(f"claim family, 50 values of t in (0, pi/2]: max ||N(gamma) - gamma|| = {worst_steady:.1e}, min TTR residual = {min_res:.3e}; at t = pi/2 residual {r_half:.4f}")
ok &= worst_steady < 1e-10 and min_res > 1e-4
U = unitary_exp(H, np.pi / 2); sg = hsqrt(g); ri = np.linalg.inv(hsqrt(ptrE(U @ np.kron(g, xi) @ U.conj().T, 3, 3)))
Nadj = lambda Y: ptrE(U.conj().T @ np.kron(Y, np.eye(3)) @ U @ np.kron(np.eye(3), xi), 3, 3)
E = lambda i: np.diag(np.eye(3)[i]).astype(complex)
R = np.array([[np.trace(sg @ Nadj(ri @ E(n) @ ri) @ sg @ E(m)).real for n in range(3)] for m in range(3)])
R_exact = np.array([[1 / 3, 17 / 26, 17 / 84], [13 / 51, 1 / 6, 13 / 28], [7 / 17, 7 / 39, 1 / 3]])
print(f"Petz population matrix at t = pi/2 (column = input) minus the exact rationals: {np.abs(R - R_exact).max():.1e}; p(1|0), p(2|1), p(0|2) = {R[1, 0]:.6f}, {R[2, 1]:.6f}, {R[0, 2]:.6f}")
ok &= np.abs(R - R_exact).max() < 1e-12
print("ALL_CHECKS_OK" if ok else "CHECK_FAIL")
sys.exit(0 if ok else 1)
