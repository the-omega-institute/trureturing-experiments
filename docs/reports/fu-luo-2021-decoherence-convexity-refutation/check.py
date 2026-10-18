#!/usr/bin/env python3
"""Check: the decoherence quantifier D(rho|Pi,U) of Fu-Luo (Entropy 23 (2021) 1594, Eq. (27)) is not convex in rho.
Part A is exact (sympy): d = 4 = 2 x 2 pointer basis, four qubit sub-environments in |0>, detectors U = (X, X, Z, Z);
the reduced state is computed literally as tr_E[Pi_U (rho (x) rho_E) Pi_U^dagger] (Eqs. (18)-(21)).
Part B is floating point (numpy): the full-rank family with all detector overlaps in (0, 1), d = 4, 6, 8, and a random
search for violations in d = 2 and d = 3. Exit 0 iff every check in Part A and Part B's family passes."""
import itertools, sys
import sympy as sp
import numpy as np

fails = []
def check(name, ok):
    print(('PASS ' if ok else 'FAIL ') + name)
    if not ok: fails.append(name)

# ---------- Part A: exact instance ----------
d = 4; p = sp.Rational(1, 4)
psi = sp.Matrix([sp.sqrt(p), sp.sqrt(1 - p)])
f = [sp.Matrix([1, 1]) / sp.sqrt(2), sp.Matrix([1, -1]) / sp.sqrt(2)]
kron = lambda A, B: sp.Matrix(sp.BlockMatrix([[A[i, j] * B for j in range(A.cols)] for i in range(A.rows)]).as_explicit())
P = psi * psi.H
rho = [kron(P, f[k] * f[k].H) for k in range(2)]
rhobar = (rho[0] + rho[1]) / 2
X = sp.Matrix([[0, 1], [1, 0]]); Z = sp.Matrix([[1, 0], [0, -1]]); I2 = sp.eye(2)
U = [X, X, Z, Z]
tau = sp.Matrix([[1, 0], [0, 0]])
def kron_list(Ms):
    out = Ms[0]
    for M in Ms[1:]: out = kron(out, M)
    return out
rhoE = kron_list([tau] * 4)
UE = [kron_list([U[i] if l == i else I2 for l in range(4)]) for i in range(4)]
Pis = [sp.Matrix(4, 4, lambda a, b: 1 if a == b == i else 0) for i in range(4)]
PiU = sum((kron(Pis[i], UE[i]) for i in range(4)), sp.zeros(64, 64))
check('Pi_U is unitary', sp.simplify(PiU * PiU.H - sp.eye(64)) == sp.zeros(64, 64))
def reduced(r):
    full = PiU * kron(r, rhoE) * PiU.H
    return sp.Matrix(4, 4, lambda a, b: sum(full[a * 16 + e, b * 16 + e] for e in range(16)))
omega = sp.Matrix(4, 4, lambda i, j: (UE[i] * rhoE * UE[j].H).trace())
print('INFO omega =', omega.tolist())
for name, r in (('rho_0', rho[0]), ('rho_1', rho[1]), ('rhobar', rhobar)):
    check(f'{name}: reduced state equals rho o Omega (Eq. (22))', sp.simplify(reduced(r) - r.multiply_elementwise(omega)) == sp.zeros(4, 4))
def is_psd(S):
    return all(sp.re(e) >= 0 and sp.simplify(sp.im(e)) == 0 for e in S.eigenvals())
def sqrt_cert(A, S, name):
    ok = sp.simplify(S * S - A) == sp.zeros(4, 4) and sp.simplify(S - S.H) == sp.zeros(4, 4) and is_psd(S)
    check(f'sqrt certificate for {name}', ok); return S
def h(S):
    return sp.nsimplify(sp.simplify(sum(S[i, i] ** 2 for i in range(4))))
Ireg = sp.diag(1, 1)
Dvals = {}
for k in range(2):
    A = rho[k]; B = reduced(A)
    SA = sqrt_cert(A, A, f'rho_{k}')   # rank-one projector
    SB = sqrt_cert(B, sp.diag(sp.sqrt(p / 2) * Ireg, sp.sqrt(1 - p) * f[k] * f[k].H), f'Phi(rho_{k})')
    Dvals[k] = sp.simplify(h(SB) - h(SA)); print(f'INFO D(rho_{k}) =', Dvals[k])
A = rhobar; B = reduced(A)
SA = sqrt_cert(A, kron(P, Ireg) / sp.sqrt(2), 'rhobar')
SB = sqrt_cert(B, sp.diag(*[sp.sqrt(B[i, i]) for i in range(4)]), 'Phi(rhobar)')
check('Phi(rhobar) is diagonal', all(B[i, j] == 0 for i in range(4) for j in range(4) if i != j))
Dbar = sp.simplify(h(SB) - h(SA)); print('INFO D(rhobar) =', Dbar)
check('D(rho_0) = D(rho_1) = 5/16', Dvals[0] == sp.Rational(5, 16) and Dvals[1] == sp.Rational(5, 16))
check('D(rhobar) = 3/8', Dbar == sp.Rational(3, 8))
check('Jensen gap D(rhobar) - (D(rho_0) + D(rho_1))/2 = 1/16 > 0', sp.simplify(Dbar - (Dvals[0] + Dvals[1]) / 2) == sp.Rational(1, 16))

# ---------- Part B: floating point ----------
def msqrt(A):
    w, V = np.linalg.eigh((A + A.conj().T) / 2); w = np.clip(w, 0, None)
    return (V * np.sqrt(w)) @ V.conj().T
def hnum(A):
    return float(np.sum(np.abs(np.diag(msqrt(A))) ** 2))
def Dnum(r, Om):
    return hnum(r * Om) - hnum(r)
def family_gap(m, pp, eps, a, b):
    dd = 2 * m
    psi = np.array([np.sqrt(pp), np.sqrt(1 - pp)]); Pp = np.outer(psi, psi.conj())
    F = np.exp(2j * np.pi * np.outer(np.arange(m), np.arange(m)) / m) / np.sqrt(m)
    rs = [(1 - eps) * np.kron(Pp, np.outer(F[:, k], F[:, k].conj())) + eps * np.eye(dd) / dd for k in range(m)]
    w = np.array([a] * m + [b] * m)
    Om = np.outer(w, w); np.fill_diagonal(Om, 1)
    rb = sum(rs) / m
    return Dnum(rb, Om) - sum(Dnum(r, Om) for r in rs) / m
for m in (2, 3, 4):
    g0 = family_gap(m, 0.25, 0.0, 0.0, 1.0)
    check(f'family m={m} (d={2*m}), degenerate detectors: gap {g0:.6f} = (1-1/m)p(1-2p) = {(1-1/m)*0.25*0.5:.6f}', abs(g0 - (1 - 1/m) * 0.125) < 1e-6)
    qd = 2.0 ** -20 / (2 * m)
    g1 = family_gap(m, 0.25, qd, qd, 1 - qd)
    check(f'family m={m} (d={2*m}), full rank, overlaps in (0,1): gap {g1:.6f} > 13/256', g1 > 13 / 256)
    g2 = family_gap(m, 0.25, 0.01, 0.1, 0.9)
    print(f'INFO family m={m}, eps=0.01, a=0.1, b=0.9: gap {g2:.6f}')
rng = np.random.default_rng(20261019)
def rand_state(dd):
    G = rng.normal(size=(dd, dd)) + 1j * rng.normal(size=(dd, dd)); r = G @ G.conj().T; return r / np.trace(r).real
for dd in (2, 3):
    worst = 0.0
    for _ in range(20000):
        om = rng.uniform(0, 1, size=dd) * np.exp(1j * rng.uniform(0, 2 * np.pi, size=dd))
        Om = np.outer(om, om.conj()); np.fill_diagonal(Om, 1)
        r1, r2 = rand_state(dd), rand_state(dd); t = rng.uniform()
        gap = Dnum(t * r1 + (1 - t) * r2, Om) - t * Dnum(r1, Om) - (1 - t) * Dnum(r2, Om)
        worst = max(worst, gap)
    print(f'INFO random search d={dd}: largest Jensen gap found {worst:.3e} (positive means a convexity violation)')
print('RESULT', 'FAIL' if fails else 'PASS', len(fails))
sys.exit(1 if fails else 0)
