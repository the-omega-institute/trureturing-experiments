#!/usr/bin/env python3
"""Alsing-Birrittella (arXiv:2404.13725v1; Foundations 5 (2025) 19), Eq. (mixed:state:summary:again): for a pure state
decomposition rho = sum_i p_i |psi_i><psi_i| with positive Hermitian amplitude matrices C_i (|psi_i> = sum C_i,nm |n,m>),
LN_a^avg(rho) <= LN_e(rho) <= LN_a(rho), with equality on the right when every C_i is real.

LN_e = log2 ||rho^Gamma||_1 (partial transpose on the second factor); LN_a(pure) = log2(1 + 2 N_a), N_a = 1/2 sum_{n, m != n}
|c_nn c_mm - c_nm c_mn|; LN_a(rho) = log2(1 + 2 N_a(rho)), N_a(rho) = 1/2 sum_{n, m != n} |rho_{nn,mm} - rho_{nm,mn}|.
Checks on 3000 random ensembles (d in {2,3,4,5}, 1-5 components, one third real) and on an exact rational example:
(1) the sandwich; (2) equality for real C_i; (3) rho^Gamma = Q F with Q = sum p_i C_i (x) C_i and F the swap;
(4) ||rho^Gamma||_1 = sum_i p_i (tr C_i)^2; (5) the exact example p = (1/2,1/2), C_1 = diag(3/5,4/5),
C_2 = [[237, 84i], [-84i, 188]]/325 gives LN values log2(119/65) < log2(7753/4225) < log2(207937/105625).
Exit 0 and RESULT PASS when every check holds within TOL.
"""
import numpy as np, sys
from fractions import Fraction as Fr
TOL = 1e-10
rng = np.random.default_rng(404)
def randC(d, real=False):
    A = rng.normal(size=(d, d)) + (0 if real else 1j * rng.normal(size=(d, d)))
    C = A @ A.conj().T; return C / np.sqrt(np.trace(C @ C).real)
def pt(R, d): return R.reshape(d, d, d, d).transpose(0, 3, 2, 1).reshape(d * d, d * d)
def tn(X): return np.linalg.svd(X, compute_uv=False).sum()
def LNa_pure(C):
    d = len(C); return np.log2(1 + 2 * sum(abs(C[n, n] * C[m, m] - C[n, m] * C[m, n]) for n in range(d) for m in range(n + 1, d)))
def LNa_mixed(R, d):
    R4 = R.reshape(d, d, d, d); return np.log2(1 + 2 * sum(abs(R4[n, n, m, m] - R4[n, m, m, n]) for n in range(d) for m in range(n + 1, d)))
e1 = e2 = e3 = e4 = 0.0
for trial in range(3000):
    d = int(rng.choice([2, 3, 4, 5])); K = int(rng.integers(1, 6)); real = trial % 3 == 0
    Cs = [randC(d, real) for _ in range(K)]; p = rng.dirichlet(np.ones(K))
    R = sum(p[k] * np.outer(Cs[k].reshape(-1), Cs[k].reshape(-1).conj()) for k in range(K))
    lo = sum(p[k] * LNa_pure(Cs[k]) for k in range(K)); ex = np.log2(tn(pt(R, d))); hi = LNa_mixed(R, d)
    e1 = max(e1, lo - ex, ex - hi)
    if real: e2 = max(e2, abs(hi - ex))
    Q = sum(p[k] * np.kron(Cs[k], Cs[k]) for k in range(K)); F = np.zeros((d * d, d * d))
    for n in range(d):
        for m in range(d): F[n * d + m, m * d + n] = 1
    e3 = max(e3, np.abs(pt(R, d) - Q @ F).max())
    e4 = max(e4, abs(tn(pt(R, d)) - sum(p[k] * np.trace(Cs[k]).real ** 2 for k in range(K))))
C1 = np.diag([0.6, 0.8]).astype(complex); C2 = np.array([[237, 84j], [-84j, 188]]) / 325
R = 0.5 * np.outer(C1.reshape(-1), C1.reshape(-1).conj()) + 0.5 * np.outer(C2.reshape(-1), C2.reshape(-1).conj())
vals = (0.5 * LNa_pure(C1) + 0.5 * LNa_pure(C2), np.log2(tn(pt(R, 2))), LNa_mixed(R, 2))
want = (np.log2(119 / 65), np.log2(7753 / 4225), np.log2(207937 / 105625))
e5 = max(abs(a - b) for a, b in zip(vals, want)); strict = vals[0] < vals[1] < vals[2]
print(f'(1) max sandwich violation = {e1:.3e}'); print(f'(2) max |LN_a(rho) - LN_e| for real C = {e2:.3e}')
print(f'(3) max |rho^Gamma - Q F| = {e3:.3e}'); print(f'(4) max |trace norm - sum p (tr C)^2| = {e4:.3e}')
print(f'(5) exact example deviation = {e5:.3e}, strict = {strict}')
ok = max(e1, e2, e3, e4, e5) < TOL and strict
print('RESULT PASS' if ok else 'RESULT FAIL'); sys.exit(0 if ok else 1)
