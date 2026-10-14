#!/usr/bin/env python3
"""Numerical checks for Rana-Parashar-Lewenstein, arXiv:1511.01854v2 (PRA 93, 012110 (2016)), Conjecture 7:
for all states rho, C_l1(rho) >= C_r(rho), with C_l1 = sum_{i != j} |rho_ij| and C_r = S2(diag rho) - S2(rho) (bits).

Proof route under test (sequential binary pinching):
  (A) Audenaert's bound for equal priors: J2(s, t) = S2((s+t)/2) - (S2(s)+S2(t))/2 <= T(s, t) = ||s - t||_tr / 2;
  (B) one pinching step: S2(P s P + Q s Q) - S2(s) <= 2 ||P s Q||_tr;
  (C) row bound: C_r(rho) <= 2 sum_k (sum_{j>k} |rho_kj|^2)^(1/2) <= C_l1(rho).
Random complex states of every rank, d = 2..8, including near-diagonal states with C_l1 < 1 (the open regime) and
rank-deficient pairs for (A). Prints the worst ratios and ALL_OK.
"""
import numpy as np
rng = np.random.default_rng(20261015)
def S2(r):
    w = np.linalg.eigvalsh((r + r.conj().T) / 2); w = w[w > 1e-15]
    return float(-(w * np.log2(w)).sum())
def rand_state(d, rank, scale=1.0):
    G = rng.normal(size=(d, rank)) + 1j * rng.normal(size=(d, rank))
    r = G @ G.conj().T; r /= np.trace(r).real
    if scale < 1.0:  # pull towards the dephased state: small coherences
        r = scale * r + (1 - scale) * np.diag(np.diag(r))
    return r
def trnorm(X): return float(np.abs(np.linalg.eigvalsh((X + X.conj().T) / 2)).sum()) if np.allclose(X, X.conj().T) else float(np.linalg.svd(X, compute_uv=False).sum())
ok = True; worst = {'A': -1e9, 'B': -1e9, 'C1': -1e9, 'C2': -1e9}
for it in range(40000):
    d = int(rng.integers(2, 9)); rank = int(rng.integers(1, d + 1))
    scale = float(rng.choice([1.0, 0.5, 0.2, 0.05, 0.01]))
    r = rand_state(d, rank, scale)
    Cl1 = float(np.abs(r).sum() - np.abs(np.diag(r)).sum())
    Cr = S2(np.diag(np.diag(r))) - S2(r)
    row = 2 * sum(np.sqrt((np.abs(r[k, k+1:]) ** 2).sum()) for k in range(d - 1))
    worst['C1'] = max(worst['C1'], Cr - row); worst['C2'] = max(worst['C2'], row - Cl1)
    # (B): a random projection
    m = int(rng.integers(1, d)); U = np.linalg.qr(rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d)))[0]
    P = U[:, :m] @ U[:, :m].conj().T; Q = np.eye(d) - P
    lhs = S2(P @ r @ P + Q @ r @ Q) - S2(r); rhs = 2 * float(np.linalg.svd(P @ r @ Q, compute_uv=False).sum())
    worst['B'] = max(worst['B'], lhs - rhs)
    # (A): a random pair, possibly rank deficient
    t = rand_state(d, int(rng.integers(1, d + 1)))
    J = S2((r + t) / 2) - (S2(r) + S2(t)) / 2; T = 0.5 * trnorm(r - t)
    worst['A'] = max(worst['A'], J - T)
for k, v in worst.items():
    print('max violation %s = %.3e' % (k, v))
    if v > 1e-9: ok = False
print('ALL_OK' if ok else 'FAIL')
