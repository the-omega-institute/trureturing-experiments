#!/usr/bin/env python3
"""Fannes-de Melo-Roga-Zyczkowski (arXiv:1104.2271, QIC 12 (2012) 472), Conjecture 1: for three density matrices
rho_i and weights p_i, chi = S(sum p_i rho_i) - sum p_i S(rho_i) <= S(C), C_ij = sqrt(p_i p_j) ||sqrt(rho_i) sqrt(rho_j)||_1.

Numerical checks of the proof route (random ensembles, d in {2,3,4,5,7}, random ranks and weights):
(1) the simultaneous rank-one decomposition rho_i = sum_a |v_ia><v_ia| built from polar factors preserves all three
    root fidelities: sum_a |<v_ia, v_ja>| = ||sqrt(rho_i) sqrt(rho_j)||_1;
(2) chi <= Xi, the Holevo quantity of the labelled ensemble;
(3) entrywise modulus of a 3x3 positive semidefinite matrix does not decrease the entropy (blocks and generic);
(4) the conjecture itself, chi <= S(C);
(5) an equality case with mixed, pairwise noncommuting states (d = 4, rho_i = |x_i><x_i| (x) I/2, x = 0, +, 1).
Exit 0 and RESULT PASS when every check holds within the floating-point tolerance TOL.
"""
import numpy as np, sys
TOL = 1e-10
rng = np.random.default_rng(20261019)
def msqrt(r):
    w, v = np.linalg.eigh((r + r.conj().T) / 2); w = np.clip(w, 0, None); return (v * np.sqrt(w)) @ v.conj().T
def polar(X):
    L, s, Rh = np.linalg.svd(X); return L @ Rh
def S(X):
    w = np.linalg.eigvalsh((X + X.conj().T) / 2); w = w[w > 1e-14]; return float(-(w * np.log(w)).sum())
def tn(X): return float(np.linalg.svd(X, compute_uv=False).sum())
def rand_state(d, r):
    A = rng.normal(size=(d, r)) + 1j * rng.normal(size=(d, r)); R = A @ A.conj().T; return R / np.trace(R).real
def unitary_eigvecs(U):
    w, V = np.linalg.eig(U)
    # U is normal: orthonormalise within each eigenvalue cluster
    out = np.zeros_like(V); used = np.zeros(len(w), bool)
    for i in range(len(w)):
        if used[i]: continue
        idx = [j for j in range(len(w)) if not used[j] and abs(w[j] - w[i]) < 1e-8]
        Q, _ = np.linalg.qr(V[:, idx]); out[:, idx] = Q; used[idx] = True
    return out
e1 = e2 = e3 = e4 = 0.0; trials = 2000
for t in range(trials):
    d = int(rng.choice([2, 3, 4, 5, 7])); rho = [rand_state(d, int(rng.integers(1, d + 1))) for _ in range(3)]
    p = rng.dirichlet(np.ones(3)); R = [msqrt(x) for x in rho]
    W = [R[0]] + [R[j] @ polar(R[0] @ R[j]).conj().T for j in (1, 2)]
    E = unitary_eigvecs(polar(W[1].conj().T @ W[2]))
    v = [[W[i] @ E[:, a] for a in range(d)] for i in range(3)]
    f = np.array([[tn(R[i] @ R[j]) for j in range(3)] for i in range(3)])
    g = np.array([[sum(abs(np.vdot(v[i][a], v[j][a])) for a in range(d)) for j in range(3)] for i in range(3)])
    e1 = max(e1, float(np.abs(f - g).max()))
    C = np.sqrt(np.outer(p, p)) * f
    chi = S(sum(p[i] * rho[i] for i in range(3))) - sum(p[i] * S(rho[i]) for i in range(3))
    sig = [sum(p[i] * np.outer(v[i][a], v[i][a].conj()) for i in range(3)) for a in range(d)]
    q = np.array([[np.vdot(v[i][a], v[i][a]).real for a in range(d)] for i in range(3)])
    Xi = sum(p[i] * q[i, a] * np.log(q[i, a]) for i in range(3) for a in range(d) if q[i, a] > 1e-15) + sum(S(s) for s in sig)
    e2 = max(e2, chi - Xi)
    for a in range(d):
        G = np.array([[np.sqrt(p[i] * p[j]) * np.vdot(v[i][a], v[j][a]) for j in range(3)] for i in range(3)])
        e3 = max(e3, S(G) - S(np.abs(G)))
    e4 = max(e4, chi - S(C))
for t in range(20000):
    A = rng.normal(size=(3, 3)) + 1j * rng.normal(size=(3, 3)); G = A @ A.conj().T; G /= np.trace(G).real
    e3 = max(e3, S(G) - S(np.abs(G)))
k0 = np.array([1, 0]); k1 = np.array([0, 1]); kp = (k0 + k1) / np.sqrt(2)
rho = [np.kron(np.outer(x, x), np.eye(2) / 2) for x in (k0, kp, k1)]; p = np.ones(3) / 3
C = np.sqrt(np.outer(p, p)) * np.array([[tn(msqrt(rho[i]) @ msqrt(rho[j])) for j in range(3)] for i in range(3)])
chi = S(sum(p[i] * rho[i] for i in range(3))) - sum(p[i] * S(rho[i]) for i in range(3))
e5 = abs(chi - S(C))
ok = e1 < TOL and e2 < TOL and e3 < TOL and e4 < TOL and e5 < TOL
print(f'(1) max |root fidelity - sum_a |<v_ia,v_ja>|| = {e1:.3e}')
print(f'(2) max chi - Xi = {e2:.3e}')
print(f'(3) max S(G) - S(|G|) = {e3:.3e}')
print(f'(4) max chi - S(C) = {e4:.3e}  (trials {trials})')
print(f'(5) equality case |chi - S(C)| = {e5:.3e}, chi = {chi:.12f}')
print('RESULT PASS' if ok else 'RESULT FAIL'); sys.exit(0 if ok else 1)
