#!/usr/bin/env python3
"""Numerical check of the Mista-Fiurasek conjecture (PRA 74, 022316 (2006); arXiv:quant-ph/0606112v1,
Eq. (2)): the optimal trade-off between the mean output fidelity F and the mean estimation fidelity G
for N copies of a Haar-random pure qudit is
    sqrt(F - 1/(N+d)) = sqrt((d-1)((N+1)/(N+d) - G)) + sqrt(N (G - N/(N+d))).

With the source's covariant reduction (Section II, Eqs. (6)-(12)), every operation is a seat chi >= 0 on
Sym^N(C^d) (x) C^d with Tr chi = D(N,d), F = Tr(chi R_F), G = Tr(chi R_G), where
    R_F = (Pi_{+,N+1})^{T_N} / D(N+1,d),
    R_G = [ (1 (x) <0|) Pi_{+,N+1} (1 (x) |0>) / D(N+1,d) ]^T (x) 1_out .
So the support function of the achievable region is h(p) = max pF + (1-p)G = D(N,d) lambda_max(p R_F + (1-p) R_G).
The check compares h(p) with the support function of the conjectured curve, for many p, and the
seat's closed form [2N+1+p(d-N)+sqrt((1-p(d-N))^2+4p^2N(d-1))]/[2(N+d)].
"""
import itertools, sys, math
import numpy as np

def sym_basis(n, d):
    """Orthonormal basis of Sym^n(C^d) inside (C^d)^{(x)n} as columns."""
    occs = [o for o in itertools.product(range(n+1), repeat=d) if sum(o) == n]
    cols = []
    for o in occs:
        v = np.zeros(d**n)
        words = set(itertools.permutations([i for i in range(d) for _ in range(o[i])]))
        for w in words:
            idx = 0
            for x in w: idx = idx*d + x
            v[idx] = 1.0
        cols.append(v / np.linalg.norm(v))
    return np.array(cols).T, occs

def run(N, d):
    BN, _ = sym_basis(N, d)          # d^N x D(N,d)
    BN1, _ = sym_basis(N+1, d)       # d^(N+1) x D(N+1,d)
    DN, DN1 = BN.shape[1], BN1.shape[1]
    P1 = BN1 @ BN1.T                 # projector on Sym^{N+1} in (C^d)^{(x)(N+1)}, real
    # restrict input to Sym^N: the input space is Sym^N (x) C^d (first N qudits symmetric, last = output)
    J = np.kron(BN, np.eye(d))       # d^(N+1) x (DN*d)
    # partial transpose on the first N qudits, then restrict
    T = P1.reshape(d**N, d, d**N, d).transpose(2, 1, 0, 3).reshape(d**(N+1), d**(N+1))
    RF = J.T @ T @ J / DN1
    e0 = np.zeros(d); e0[0] = 1
    K = np.kron(np.eye(d**N), e0.reshape(d, 1))
    RGin = BN.T @ (K.T @ P1 @ K) @ BN / DN1
    RG = np.kron(RGin.T, np.eye(d))
    worst = 0.0
    for p in np.linspace(0, 1, 41):
        h = DN * np.linalg.eigvalsh(p*RF + (1-p)*RG).max()
        closed = (2*N + 1 + p*(d - N) + math.sqrt((1 - p*(d - N))**2 + 4*p*p*N*(d - 1))) / (2*(N + d))
        # support of the conjectured curve
        Gs = np.linspace(N/(N+d-1), (N+1)/(N+d), 4001)
        Fs = 1/(N+d) + (np.sqrt(np.clip((d-1)*((N+1)/(N+d) - Gs), 0, None)) + np.sqrt(np.clip(N*(Gs - N/(N+d)), 0, None)))**2
        hc = np.max(p*Fs + (1-p)*Gs)
        worst = max(worst, abs(h - closed), abs(h - hc) if p > 0 else 0)
    return worst, DN

ok = True
for d in (2, 3, 4):
    for N in (1, 2, 3, 4):
        if d**(N+1) > 5000: continue
        w, DN = run(N, d)
        line = w < 2e-6
        ok &= line
        print(f"d={d} N={N} dim Sym^N={DN}: max |support - closed form or curve support| = {w:.2e} {'OK' if line else 'FAIL'}")
print("ALL_OK" if ok else "SOME_FAIL")
sys.exit(0 if ok else 1)
