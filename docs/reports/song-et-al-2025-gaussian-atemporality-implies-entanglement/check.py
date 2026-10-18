#!/usr/bin/env python3
"""Check: for physical two-mode Gaussian covariance matrices, PPT implies forward and reverse temporal compatibility
(Song et al., arXiv:2508.11804v2, Eqs. (5), (7), (8)); hence Gaussian atemporality implies entanglement.
Part A (sympy, symbolic): the congruence identities K V K^T = N, K((-Ω)⊕Ω)K^T = Ω - TΩT^T and the reverse ones.
Part B (numpy): random physical standard-form covariance matrices; every PPT one is compatible in both directions,
every atemporal one found is NPT; the source's beam-splitter family. Exit 0 iff all checks pass."""
import sys
import sympy as sp
import numpy as np
fails = []
def check(name, ok):
    print(('PASS ' if ok else 'FAIL ') + name)
    if not ok: fails.append(name)
# ---------- Part A ----------
a1, a2, b1, b2 = sp.symbols('a1 a2 b1 b2', positive=True)
c11, c12, c21, c22 = sp.symbols('c11 c12 c21 c22', real=True)
A = sp.diag(a1, a2); B = sp.diag(b1, b2); C = sp.Matrix([[c11, c12], [c21, c22]])
Om = sp.Matrix([[0, 1], [-1, 0]])
V = sp.BlockMatrix([[A, C], [C.T, B]]).as_explicit()
T = C.T * A.inv(); N = B - C.T * A.inv() * C
K = sp.Matrix.hstack(-T, sp.eye(2))
mOmOm = sp.diag(-Om, Om)
check('K V K^T = N', sp.simplify(K * V * K.T - N) == sp.zeros(2, 2))
check('K((-Ω)⊕Ω)K^T = Ω - TΩT^T', sp.simplify(K * mOmOm * K.T - (Om - T * Om * T.T)) == sp.zeros(2, 2))
check('A T^T = C and T A T^T + N = B (Eq. (6))', sp.simplify(A * T.T - C) == sp.zeros(2, 2) and sp.simplify(T * A * T.T + N - B) == sp.zeros(2, 2))
S = C * B.inv(); M = A - C * B.inv() * C.T
L = sp.Matrix.hstack(sp.eye(2), -S)
OmmOm = sp.diag(Om, -Om)
check('L V L^T = M', sp.simplify(L * V * L.T - M) == sp.zeros(2, 2))
check('L(Ω⊕(-Ω))L^T = Ω - SΩS^T', sp.simplify(L * OmmOm * L.T - (Om - S * Om * S.T)) == sp.zeros(2, 2))
Lam = sp.diag(1, -1, 1, 1)
check('Λ(V + iΩ⊕Ω)Λ = ΛVΛ + i((-Ω)⊕Ω) (partial transpose flips p_A)', sp.simplify(Lam * sp.diag(Om, Om) * Lam - mOmOm) == sp.zeros(4, 4))
# ---------- Part B ----------
O = np.array([[0, 1], [-1, 0]])
def psd(H, tol=1e-9): return np.linalg.eigvalsh((H + H.conj().T) / 2).min() >= -tol
def parts(V):
    A, C, B = V[:2, :2], V[:2, 2:], V[2:, 2:]
    T = C.T @ np.linalg.inv(A); N = B - C.T @ np.linalg.inv(A) @ C
    S = C @ np.linalg.inv(B); M = A - C @ np.linalg.inv(B) @ C.T
    Xf = N + 1j * O - 1j * T @ O @ T.T; Xr = M + 1j * O - 1j * S @ O @ S.T
    return Xf, Xr
def physical(V): return psd(V + 1j * np.kron(np.eye(2), O))
def ppt(V):
    Lam = np.diag([1, -1, 1, 1]); return psd(Lam @ V @ Lam + 1j * np.kron(np.eye(2), O))
rng = np.random.default_rng(20261019)
def symplectic_random():
    # product of random two-mode squeezers, beam splitters and phase rotations
    def rot(t): return np.array([[np.cos(t), np.sin(t)], [-np.sin(t), np.cos(t)]])
    Sm = np.eye(4)
    for _ in range(4):
        r = rng.normal(scale=0.7, size=2); sq = np.diag([np.exp(r[0]), np.exp(-r[0]), np.exp(r[1]), np.exp(-r[1])])
        th = rng.uniform(0, 2 * np.pi); bs = np.block([[np.cos(th) * np.eye(2), np.sin(th) * np.eye(2)], [-np.sin(th) * np.eye(2), np.cos(th) * np.eye(2)]])
        ph = np.block([[rot(rng.uniform(0, 6.3)), np.zeros((2, 2))], [np.zeros((2, 2)), rot(rng.uniform(0, 6.3))]])
        Sm = ph @ bs @ sq @ Sm
    return Sm
def to_standard(V):
    # local symplectic diagonalisation of each block (Williamson-free: rotate then squeeze each mode to a diagonal block)
    out = V.copy()
    for k in (0, 2):
        blk = out[k:k + 2, k:k + 2]; w, U = np.linalg.eigh(blk)
        R = U.T if np.linalg.det(U) > 0 else (U @ np.diag([1, -1])).T
        F = np.eye(4); F[k:k + 2, k:k + 2] = R; out = F @ out @ F.T
    return out
counts = dict(total=0, ppt=0, ppt_both=0, atemporal=0, atemporal_npt=0, oneway=0, oneway_npt=0)
for _ in range(20000):
    Sm = symplectic_random(); nu = 1 + rng.exponential(0.5, size=2)
    V = Sm @ np.diag([nu[0], nu[0], nu[1], nu[1]]) @ Sm.T
    V = to_standard(V)
    if not physical(V): continue
    counts['total'] += 1
    Xf, Xr = parts(V); f, r = psd(Xf), psd(Xr)
    if ppt(V):
        counts['ppt'] += 1; counts['ppt_both'] += int(f and r)
    if (not f) and (not r):
        counts['atemporal'] += 1; counts['atemporal_npt'] += int(not ppt(V))
    if (not f) or (not r):
        counts['oneway'] += 1; counts['oneway_npt'] += int(not ppt(V))
print('INFO', counts)
check('every PPT sample is compatible in both directions', counts['ppt'] == counts['ppt_both'] and counts['ppt'] > 0)
check('every atemporal sample is NPT (entangled)', counts['atemporal'] == counts['atemporal_npt'] and counts['atemporal'] > 0)
check('every sample incompatible in at least one direction is NPT', counts['oneway'] == counts['oneway_npt'])
# The source's beam-splitter family: atemporal iff (1-u)/u^2 > (v-1)/v^2 (Appendix B), entangled iff min(u,v) < 1.
bad = 0; seen = 0
for u in np.linspace(0.2, 0.99, 40):
    for v in np.linspace(1.01, 6, 40):
        if u * v < 1: continue
        V = 0.5 * np.array([[u + v, 0, -(u - v), 0], [0, u + v, 0, u - v], [-(u - v), 0, u + v, 0], [0, u - v, 0, u + v]])
        Xf, Xr = parts(V)
        if not psd(Xf) and not psd(Xr):
            seen += 1; bad += int(ppt(V))
print(f'INFO beam-splitter family: {seen} atemporal grid points')
check('beam-splitter family: atemporal points are NPT', bad == 0 and seen > 0)
print('RESULT', 'FAIL' if fails else 'PASS', len(fails))
sys.exit(1 if fails else 0)
