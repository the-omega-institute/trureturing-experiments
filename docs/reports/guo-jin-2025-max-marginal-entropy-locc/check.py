#!/usr/bin/env python3
"""Checks for Y. Guo and Z.-X. Jin, "Measure of entanglement and the monogamy relation: a topical review",
arXiv:2512.21992v1, Section 9.4: the conjecture that E'^(n)(psi) = max_i h(rho_{A_i}) (Eq. (189)) does not increase on
average under LOCC, for h = S (von Neumann entropy, -Tr rho log2 rho) and the other reduced functions of Table 16.

Fidelity: E'_f is "proper" (the definition before Eq. (157)): E'_f(GHZ_n) = 1 > E'_f(W_n) = H2(1/n), n = 3..8.
Claim: psi_d = (1/sqrt(2d)) sum_{j=1}^d (|j,j,0,0> + |0,0,j,j>) on (C^{d+1})^{(x)4}; Alice's two-outcome projective
measurement P0 = |0><0|, P1 = 1 - P0 raises the average of max_i S(rho_i) from 1 + (log2 d)/2 to log2 d (d = 8: 5/2 -> 3).
Also for h = tau = 2(1 - Tr rho^2) and h = 1 - ||rho||. Observation: for n = 3, a random search over qubit and qutrit
states with random two-outcome local projective measurements finds no increase. Exits 1 on failure.
"""
import math, sys
import numpy as np
ok = True
def ent(p):
    w = np.linalg.eigvalsh((p + p.conj().T) / 2); w = w[w > 1e-14]; return float(-(w * np.log2(w)).sum())
def tau(p): return float(2 * (1 - np.trace(p @ p).real))
def opn(p): return float(1 - np.linalg.eigvalsh((p + p.conj().T) / 2).max())
def marg(psi, dims, i):
    t = psi.reshape(dims); t = np.moveaxis(t, i, 0).reshape(dims[i], -1); return t @ t.conj().T
def Emax(psi, dims, h=ent): return max(h(marg(psi, dims, i)) for i in range(len(dims)))
def H2(x): return -x * math.log2(x) - (1 - x) * math.log2(1 - x)
fid = True
for n in range(3, 9):
    dims = (2,) * n
    ghz = np.zeros(2 ** n); ghz[0] = ghz[-1] = 1 / math.sqrt(2)
    w = np.zeros(2 ** n)
    for k in range(n): w[1 << k] = 1 / math.sqrt(n)
    g, ww = Emax(ghz, dims), Emax(w, dims)
    fid &= abs(g - 1) < 1e-12 and abs(ww - H2(1 / n)) < 1e-12 and g > ww
print(f"fidelity: E'_f(GHZ_n) = 1 > E'_f(W_n) = H2(1/n) for n = 3..8: {fid}")
ok &= fid
def family(d):
    D = d + 1; dims = (D,) * 4
    psi = np.zeros(D ** 4)
    idx = lambda a, b, c, e: ((a * D + b) * D + c) * D + e
    for j in range(1, d + 1):
        psi[idx(j, j, 0, 0)] = psi[idx(0, 0, j, j)] = 1 / math.sqrt(2 * d)
    P0 = np.zeros((D, D)); P0[0, 0] = 1; P1 = np.eye(D) - P0
    outs = []
    for P in (P0, P1):
        v = np.kron(P, np.eye(D ** 3)) @ psi; p = float(v @ v); outs.append((p, v / math.sqrt(p)))
    return psi, dims, outs
for d in (8, 5, 6, 7, 9, 12, 16):
    psi, dims, outs = family(d)
    init = Emax(psi, dims); avg = sum(p * Emax(v, dims) for p, v in outs)
    good = abs(np.linalg.norm(psi) - 1) < 1e-12 and abs(init - (1 + math.log2(d) / 2)) < 1e-9 and abs(avg - math.log2(d)) < 1e-9 and avg > init
    margs = [np.round(np.diag(marg(psi, dims, i)).real, 12) for i in range(4)]
    offdiag = max(np.abs(marg(psi, dims, i) - np.diag(np.diag(marg(psi, dims, i)))).max() for i in range(4))
    t_init = Emax(psi, dims, tau); t_avg = sum(p * Emax(v, dims, tau) for p, v in outs)
    o_init = Emax(psi, dims, opn); o_avg = sum(p * Emax(v, dims, opn) for p, v in outs)
    print(f"d = {d:2d}: probabilities {[round(p, 12) for p, _ in outs]}; marginals diagonal (max off-diagonal {offdiag:.1e}); "
          f"E'_f {init:.6f} -> average {avg:.6f} (increase {avg - init:.6f}); tau {t_init:.6f} -> {t_avg:.6f}; 1-||rho|| {o_init:.6f} -> {o_avg:.6f}")
    good &= offdiag < 1e-14 and t_avg > t_init and o_avg > o_init
    ok &= good
rng = np.random.default_rng(20261024)
best = -1.0
for d in (2, 3):
    dims = (d, d, d)
    for trial in range(3000):
        psi = rng.normal(size=d ** 3) + 1j * rng.normal(size=d ** 3)
        if trial % 2: psi = psi * (rng.random(d ** 3) < 0.4)
        if np.linalg.norm(psi) < 1e-9: continue
        psi /= np.linalg.norm(psi)
        Z = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d)); Q, _ = np.linalg.qr(Z)
        for k in range(1, d):
            P0 = Q[:, :k] @ Q[:, :k].conj().T
            avg = 0.0
            for P in (P0, np.eye(d) - P0):
                v = np.kron(P, np.eye(d * d)) @ psi; p = np.vdot(v, v).real
                if p > 1e-12: avg += p * Emax(v / math.sqrt(p), dims)
            best = max(best, avg - Emax(psi, dims))
print(f"observation n = 3: largest average increase of E'_f over 6000 random qubit/qutrit states and measurements: {best:.2e}")
print("ALL_CHECKS_OK" if ok else "CHECK_FAIL")
sys.exit(0 if ok else 1)
