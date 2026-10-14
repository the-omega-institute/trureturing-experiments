#!/usr/bin/env python3
"""Check of the proposed proof of the Yu-Duan-Xu equality (arXiv:1201.1172v1, IEEE TIT 63 (2017)
1299-1310, Sec. V after Theorem 3: "in Eq. (8) ... the second one should be an equality. However, this
is still an unsettled issue"): D(Phi_C, Conv(U_d)) = D(Phi_C, Lambda_d) for every Schur channel Phi_C.

The superchannel Theta(T) = R o (T (x) id_d) o Ad_V, with V|i> = |i,i> and R having Kraus operators
K0 = sum_i |i><i,i| and K_{k,i} = |i><k,i| (k != i).

Checks for d = 1..5 (float64, random instances):
 1. V^dagger V = I_d and sum K^dagger K = I_{d^2} (R and Ad_V are CPTP).
 2. Theta(T)(E_ij) = <i|T(E_ij)|j> E_ij for i != j and Theta(T)(E_ii) = Tr T(E_ii) E_ii, for random
    linear T; hence Theta(Phi_C) = Phi_C for random correlation matrices C.
 3. Theta(Ad_U) = Phi_B with B = a a^dagger + diag(1 - |a_i|^2), a_i = U_ii, and
    Phi_B = 2^-d sum_eps Ad_{D_eps} with D_eps = diag(z_i^{eps_i}), z_i^pm = q_i (r_i +- i sqrt(1 - r_i^2)).
"""
import itertools, sys
import numpy as np
rng = np.random.default_rng(20261014)

def apply_kraus(Ks, X):
    return sum(K @ X @ K.conj().T for K in Ks)

def theta(T, d):
    """Return the linear map Theta(T) as a function on d x d matrices."""
    V = np.zeros((d*d, d), dtype=complex)
    for i in range(d): V[i*d + i, i] = 1
    Ks = []
    K0 = np.zeros((d, d*d), dtype=complex)
    for i in range(d): K0[i, i*d + i] = 1
    Ks.append(K0)
    for k in range(d):
        for i in range(d):
            if k != i:
                K = np.zeros((d, d*d), dtype=complex); K[i, k*d + i] = 1; Ks.append(K)
    def out(X):
        Y = V @ X @ V.conj().T                      # Ad_V
        Z = np.zeros((d*d, d*d), dtype=complex)     # (T (x) id_d)
        for p in range(d):
            for q in range(d):
                blk = Y[np.ix_([p*d + s for s in range(d)], [q*d + s for s in range(d)])]
                # Y = sum_{p,q} E_pq (x) blk_pq on H (x) H with first factor the T-input
                Epq = np.zeros((d, d), dtype=complex); Epq[p, q] = 1
                Z += np.kron(T(Epq), blk)
        return apply_kraus(Ks, Z)
    return out, V, Ks

def rand_unitary(d):
    Z = rng.normal(size=(d, d)) + 1j*rng.normal(size=(d, d))
    Q, R = np.linalg.qr(Z); return Q @ np.diag(np.diag(R) / abs(np.diag(R)))

ok = True; worst = 0.0
for d in range(1, 6):
    # 1. CPTP
    th, V, Ks = theta(lambda X: X, d)
    e1 = max(np.abs(V.conj().T @ V - np.eye(d)).max(), np.abs(sum(K.conj().T @ K for K in Ks) - np.eye(d*d)).max())
    # 2. action on matrix units for a random linear map T(X) = sum_m A_m X B_m
    A = [rng.normal(size=(d, d)) + 1j*rng.normal(size=(d, d)) for _ in range(3)]
    Bm = [rng.normal(size=(d, d)) + 1j*rng.normal(size=(d, d)) for _ in range(3)]
    T = lambda X: sum(a @ X @ b for a, b in zip(A, Bm))
    thT, _, _ = theta(T, d); e2 = 0.0
    for i in range(d):
        for j in range(d):
            E = np.zeros((d, d), dtype=complex); E[i, j] = 1
            want = (T(E)[i, j] if i != j else np.trace(T(E))) * E
            e2 = max(e2, np.abs(thT(E) - want).max())
    # Schur channel fixed
    G = rng.normal(size=(d, d)) + 1j*rng.normal(size=(d, d)); C = G @ G.conj().T
    s = np.sqrt(np.real(np.diag(C))); C = C / np.outer(s, s)
    PhiC = lambda X: C * X
    thC, _, _ = theta(PhiC, d); X = rng.normal(size=(d, d)) + 1j*rng.normal(size=(d, d))
    e3 = np.abs(thC(X) - PhiC(X)).max()
    # 3. unitary channel -> mixture of diagonal unitaries
    U = rand_unitary(d); AdU = lambda X: U @ X @ U.conj().T
    thU, _, _ = theta(AdU, d)
    a = np.diag(U); r = np.abs(a); q = np.where(r > 0, a / np.where(r > 0, r, 1), 1)
    zp = q * (r + 1j*np.sqrt(np.clip(1 - r**2, 0, None))); zm = q * (r - 1j*np.sqrt(np.clip(1 - r**2, 0, None)))
    mix = lambda X: sum(np.diag([zp[i] if e[i] else zm[i] for i in range(d)]) @ X @ np.diag([zp[i] if e[i] else zm[i] for i in range(d)]).conj().T
                        for e in itertools.product([0, 1], repeat=d)) / 2**d
    Bmat = np.outer(a, a.conj()) + np.diag(1 - r**2)
    e4 = max(np.abs(thU(X) - mix(X)).max(), np.abs(thU(X) - Bmat * X).max())
    line = max(e1, e2, e3, e4) < 1e-12
    ok &= line; worst = max(worst, e1, e2, e3, e4)
    print(f"d={d}: CPTP residual {e1:.1e}, matrix-unit action {e2:.1e}, Schur fixed {e3:.1e}, "
          f"unitary -> diagonal mixture {e4:.1e} {'OK' if line else 'FAIL'}")
print(f"max residual {worst:.1e}")
print("ALL_OK" if ok else "SOME_FAIL")
sys.exit(0 if ok else 1)
