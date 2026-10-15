"""Exact and numerical checks for the reweighting question of Drotos-Pal-Taoutioui-Vertesi, arXiv:2406.08661 (NJP 26,
063012 (2024)), Section V: for Bloch vectors with no strictly positive null combination, can positive r_x always make
every degenerate (fixed-output) measurement inferior?  Witness W = sum_xy w_xy [mu_y + (1-|mu_y|) a_x . v_y] (Eq. (11)),
w_xy = r_x mu_xy with mu_xy the components of m_x in an eigenframe of S = sum_x r_x m_x m_x^T.  Prints ALL_OK."""
from fractions import Fraction as Fr
import numpy as np, itertools
ok = True
m = [[Fr(s*2, 7), Fr(t*3, 7), Fr(6, 7)] for s, t in [(1, 1), (1, -1), (-1, 1), (-1, -1)]]
ok &= all(sum(c*c for c in v) == 1 for v in m)                      # unit Bloch vectors
ok &= all(sum(a*b for a, b in zip(u, v)) > 0 for u in m for v in m)  # acute ensemble
# no strictly positive null combination: third coordinates are all positive
ok &= all(v[2] > 0 for v in m)
# r = 1: S = diag(16, 36, 144)/49, frame = standard axes, w = m
S = [[sum(v[i]*v[j] for v in m) for j in range(3)] for i in range(3)]
ok &= S == [[Fr(16, 49), 0, 0], [0, Fr(36, 49), 0], [0, 0, Fr(144, 49)]]
R = 4
genuine_target = sum(m[x][y]*m[x][y] for x in range(4) for y in range(3))
with_degenerate = sum(m[x][y]*m[x][y] for x in range(4) for y in range(2)) + sum(m[x][2] for x in range(4))
print("target genuine score", genuine_target, " with Z -> I:", with_degenerate)
ok &= genuine_target == R and with_degenerate == Fr(220, 49)
# numerical sanity of the global genuine optimum R: random search over unit b_y and optimal preparations
rng = np.random.default_rng(1)
M = np.array([[float(c) for c in v] for v in m])
best = 0.0
for _ in range(20000):
    B = rng.normal(size=(3, 3)); B /= np.linalg.norm(B, axis=1, keepdims=True)
    h = M @ B                       # h_x = sum_y w_xy b_y with w = M
    best = max(best, np.linalg.norm(h, axis=1).sum())
print("random genuine optimum estimate", best, "<= 4:", best <= 4 + 1e-9)
ok &= best <= 4 + 1e-9
# all positive weights and every eigenframe: Perron axis gives a degenerate improvement (random r)
for _ in range(2000):
    r = rng.uniform(0.01, 10, size=4)
    Sm = (M.T * r) @ M
    lam, V = np.linalg.eigh(Sm)
    q = M @ V                                   # q_xy
    gain = max(abs((r * q[:, y]).sum()) - (r * q[:, y]**2).sum() for y in range(3))
    ok &= gain > 0
print("ALL_OK" if ok else "FAIL")
