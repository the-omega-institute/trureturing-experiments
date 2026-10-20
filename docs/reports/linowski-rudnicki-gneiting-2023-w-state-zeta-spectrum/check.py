#!/usr/bin/env python3
"""Linowski-Rudnicki-Gneiting, arXiv:2212.12245v2 (PRA 107, 042218 (2023)), Example 4: zeta = (1/N) sum_j |zeta_j><zeta_j|,
|zeta_j> = (N-1)^(-1/2) sum_{k != j} |e_{j,k}> on (C^2)^{(x)N}. Conjectures: lambda_max(zeta) = 2/N (after Eq. (52)); and the
spectrum is 2/N (x1), (N-2)/(N(N-1)) (x(N-1)), 0 (x(2^N - N)) (note [28]).
Exact part (Fractions): F e_j = sum_{k != j} e_{j,k}; checks F^T F = (N-2) I + J and zeta = F F^T / (N(N-1)) entrywise for
N = 2..12. Floating part: eigenvalues of zeta for N = 2..10 against the conjectured spectrum. usage: check.py"""
from fractions import Fraction as Q
import itertools, sys
import numpy as np

def pairs(N): return [frozenset(p) for p in itertools.combinations(range(N), 2)]

ok = True
for N in range(2, 13):
    P = pairs(N); col = {j: [p for p in P if j in p] for j in range(N)}
    gram_ok = all((len(set(col[j]) & set(col[k])) == (N - 1 if j == k else 1)) for j in range(N) for k in range(N))
    # zeta restricted to the weight-2 sector (it vanishes elsewhere): zeta[p][q] = (1/(N(N-1))) #{j : j in p and j in q}
    zeta_ok = True
    for p in P:
        for q in P:
            direct = sum(Q(1, N) * Q(1, N - 1) for j in range(N) if j in p and j in q)
            via_F = Q(len(p & q), N * (N - 1))
            zeta_ok &= (direct == via_F)
    good = gram_ok and zeta_ok; ok &= good
    print(f"N={N}: F^T F = (N-2)I + J {gram_ok}; zeta = F F^T/(N(N-1)) on the weight-2 sector {zeta_ok}")
for N in range(2, 11):
    n = 2 ** N; Z = np.zeros((n, n))
    idx = lambda s: sum(1 << i for i in s)
    for j in range(N):
        v = np.zeros(n)
        for k in range(N):
            if k != j: v[idx({j, k})] = 1 / np.sqrt(N - 1)
        Z += np.outer(v, v) / N
    w = np.sort(np.linalg.eigvalsh(Z))[::-1]
    want = np.sort(np.array([2 / N] + [(N - 2) / (N * (N - 1))] * (N - 1) + [0.0] * (n - N)))[::-1]
    err = float(np.max(np.abs(w - want))); good = err < 1e-10; ok &= good
    print(f"N={N}: lambda_max = {w[0]:.12f} (2/N = {2 / N:.12f}); max |spectrum - conjectured| = {err:.2e}")
print("RESULT PASS" if ok else "RESULT FAIL"); print("EXIT=0" if ok else "EXIT=1"); sys.exit(0 if ok else 1)
