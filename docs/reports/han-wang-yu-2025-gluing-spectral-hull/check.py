#!/usr/bin/env python3
"""Han-Wang-Yu (Mathematics 13 (2025) 2080), Conjecture 1: the eigenvalues of a theta-gluing (Definition 8) lie in the
convex hull of the eigenvalues of the glued chains.

Exact (fractions): the seed B = (1/54)[[43,9,2],[2,36,16],[9,9,36]] (root = state 0) is column stochastic with
charpoly (x-1)(x-17/27)(x-1/2); for k = 2..4 copies and several theta (uniform and nonuniform), the glued matrix is
column stochastic and has the eigenvector (0, x, -x, 0, ...) with eigenvalue 4/9 < 1/2; charpoly of the glue equals
chi_B * chi_D^(k-1) (numerically). Also the family B_{r,delta} for (r, delta) in a grid. Exit 0 and RESULT PASS.
"""
import sys
from fractions import Fraction as F
import numpy as np

def glue(chains, theta):
    # chains: list of square matrices (lists of Fractions), root index 0
    idx = [(None, None)] + [(i, v) for i, A in enumerate(chains) for v in range(1, len(A))]
    N = len(idx); G = [[F(0)] * N for _ in range(N)]
    for a, (i, v) in enumerate(idx):
        for b, (j, w) in enumerate(idx):
            if i is None and j is None: G[a][b] = sum(theta[t] * chains[t][0][0] for t in range(len(chains)))
            elif j is None: G[a][b] = theta[i] * chains[i][v][0]
            elif i is None: G[a][b] = chains[j][0][w]
            elif i == j: G[a][b] = chains[i][v][w]
    return G, idx

def colstoch(M): return all(x >= 0 for r in M for x in r) and all(sum(M[v][w] for v in range(len(M))) == 1 for w in range(len(M)))
ok = True
B = [[F(43, 54), F(9, 54), F(2, 54)], [F(2, 54), F(36, 54), F(16, 54)], [F(9, 54), F(9, 54), F(36, 54)]]
ok &= colstoch(B)
Bn = np.array(B, dtype=float); ev = sorted(np.linalg.eigvals(Bn).real)
ok &= np.allclose(ev, [0.5, 17 / 27, 1.0])
for lam in [F(1), F(17, 27), F(1, 2)]:
    M = [[(F(1) if r == c else F(0)) * lam - B[r][c] for c in range(3)] for r in range(3)]
    d = (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
         + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))
    ok &= (d == 0)
for k in [2, 3, 4]:
    for theta in ([F(1, k)] * k, [F(1, 3), F(2, 3)] if k == 2 else [F(1, 2)] + [F(1, 2 * (k - 1))] * (k - 1)):
        G, idx = glue([B] * k, theta); ok &= colstoch(G)
        x = [F(0)] * len(idx); x[1], x[2] = F(-4), F(3); x[3], x[4] = F(4), F(-3)
        Gx = [sum(G[a][b] * x[b] for b in range(len(idx))) for a in range(len(idx))]
        ok &= all(Gx[a] == F(4, 9) * x[a] for a in range(len(idx)))
        evG = np.linalg.eigvals(np.array(G, dtype=float))
        ok &= bool(np.min(np.abs(evG - 4 / 9)) < 1e-9)
        print(f"k={k} theta={[str(t) for t in theta]} stochastic glue, eigenvalue 4/9 present; hull of inputs = [1/2, 1]")
# family B_{r,delta}
for r in [F(1, 9), F(1, 20), F(1, 100)]:
    h = 2 * r + 3 * r * r
    Br = [[F(1, 3) + h, F(1, 3), F(1, 3) - h], [F(1, 3) - h, F(1, 3), F(1, 3) + h], [F(1, 3)] * 3]
    for delta in [F(1, 2), F(3, 4), F(9, 10)]:
        Bd = [[delta * (r_ == c) + (1 - delta) * Br[r_][c] for c in range(3)] for r_ in range(3)]
        ok &= colstoch(Bd) and all(x > 0 for row in Bd for x in row)
        G, idx = glue([Bd, Bd], [F(1, 2), F(1, 2)])
        mu = delta - (1 - delta) * r
        ev = np.linalg.eigvals(np.array(G, dtype=float))
        ok &= bool(np.min(np.abs(ev - float(mu))) < 1e-9) and mu < delta
print("RESULT", "PASS" if ok else "FAIL"); sys.exit(0 if ok else 1)
