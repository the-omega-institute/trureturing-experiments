#!/usr/bin/env python3
"""Checks for Cobucci-Tavakoli, arXiv:2402.06234v4 (Sci. Adv. 10, eadq4467 (2024)), Supplementary V, Eqs. (V1)-(V5):
'We have shown this explicitly for odd d <= 17 and conjecture it to hold for odd d in general.'
W_d = 3 sum_l |lll><lll| + sum_{(j1,j2,j3) in L} sum_{l1+l2+l3=0} (x) |e_lk^(jk)><e_lk^(jk)|, e_l^(j) = d^-1/2 sum_m w^{m(l+jm)}|m>.
1. Dense numerics (numpy eigvalsh), odd d = 3, 5, 7, 9, 11, 13, 15 (composite 9, 15 included): the GHZ vector is an
   eigenvector with eigenvalue 12, eigenvalue 12 is simple, and every other eigenvalue lies in [0, 3] (so (V5) holds);
   the sorted dense spectrum equals the sorted closed-form list lambda(x,y,s) below, entry by entry, to 1e-8.
2. The closed-form spectrum lambda(x,y,s) = 3[x=y=0] + #{(a,b) in {-1,0,1}^2 : s = 2(ax+by) mod d} has max 3 off (0,0,0)
   for every odd d <= 201 (exact integer counting).
Prints ALL_OK.
"""
import numpy as np, itertools, sys
ok = True
def W(d):
    w = np.exp(2j * np.pi / d); m = np.arange(d)
    E = {j: np.array([[w ** ((mm * (l + j * mm)) % d) for mm in m] for l in m]) / np.sqrt(d) for j in range(d)}  # rows: e_l
    L = [(0,0,0),(1,1,d-2),(d-1,d-1,2),(0,1,d-1),(1,0,d-1),(0,d-1,1),(1,d-1,0),(d-1,0,1),(d-1,1,0)]
    D3 = d ** 3; M = np.zeros((D3, D3), dtype=complex)
    for l in range(d):
        idx = l * d * d + l * d + l; M[idx, idx] += 3
    for (j1, j2, j3) in L:
        for l1 in range(d):
            for l2 in range(d):
                l3 = (-l1 - l2) % d
                v = np.kron(np.kron(E[j1][l1], E[j2][l2]), E[j3][l3])
                M += np.outer(v, v.conj())
    return M
for d in (3, 5, 7, 9, 11, 13, 15):
    M = W(d)
    g = np.zeros(d ** 3, dtype=complex)
    for t in range(d): g[t * d * d + t * d + t] = 1 / np.sqrt(d)
    resid = np.linalg.norm(M @ g - 12 * g)
    ev = np.sort(np.linalg.eigvalsh((M + M.conj().T) / 2))
    top, second, low = ev[-1], ev[-2], ev[0]
    print('d=%2d  |W g - 12 g|=%.2e  top=%.6f  second=%.6f  min=%.2e' % (d, resid, top, second, low), flush=True)
    if resid > 1e-8 or abs(top - 12) > 1e-8 or second > 3 + 1e-8 or low < -1e-8: ok = False
    closed = sorted(3 * (x == 0 and y == 0) + sum(1 for a in (-1, 0, 1) for b in (-1, 0, 1) if (2 * (a * x + b * y) - s) % d == 0)
                    for x in range(d) for y in range(d) for s in range(d))
    gap = float(np.max(np.abs(ev - np.array(closed, dtype=float))))
    print('      closed-form spectrum vs dense: max |difference| = %.2e' % gap, flush=True)
    if gap > 1e-8: ok = False
for d in range(3, 202, 2):
    worst = 0
    for x in range(d):
        for y in range(d):
            if x == 0 and y == 0: continue
            cnt = {}
            for a in (-1, 0, 1):
                for b in (-1, 0, 1):
                    s = (2 * (a * x + b * y)) % d; cnt[s] = cnt.get(s, 0) + 1
            worst = max(worst, max(cnt.values()))
    if worst > 3: print('COUNT_FAIL', d, worst); ok = False
print('closed-form count: max eigenvalue off (0,0,0) <= 3 for all odd d <= 201')
print('ALL_OK' if ok else 'FAIL')
