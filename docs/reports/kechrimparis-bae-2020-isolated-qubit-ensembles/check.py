#!/usr/bin/env python3
"""Kechrimparis-Bae, New J. Phys. 22 (2020) 123024 (arXiv:2006.16366v2), Conclusions: "We expect that no isolated
ensembles exist". An ensemble is isolated when the identity is the only CPTP qubit channel preserving one of its
optimal minimum-error POVMs. Checks (not a proof) of the four-state family with G > 18 (rational instance G = 19):
  A. exact rational: states and priors are valid, M is a POVM, the dual certificate K = G I/(2L) gives
     K - q_i rho_i = r_i P_i >= 0 with P_i M_i = 0, so M is optimal with success G/L;
  B. Holevo-Yuen-Kennedy-Lax test of M on the output ensemble {q_i, N(rho_i)} fails for every sampled non-identity
     CPTP channel (random Kraus channels near and far from the identity, unitaries, depolarizing, dephasing,
     amplitude damping) and passes for the identity;
  C. the Choi matrix of the affine form forced by optimality (Lemma 2) and the two fixed quadratic forms of Lemma 3:
     v C v + w C w = -33 eps / 5, symbolically.
Exit 0 and final line `bad= 0` when every check agrees."""
import sys
from fractions import Fraction as Fr
import numpy as np
import sympy as sp

bad = 0
X = np.array([[0, 1], [1, 0]], dtype=complex); Y = np.array([[0, -1j], [1j, 0]]); Z = np.array([[1, 0], [0, -1]], dtype=complex); I2 = np.eye(2)
s = [np.array([3/5, 0, 4/5]), np.array([-3/5, 0, 4/5]), np.array([0, 3/5, -4/5]), np.array([0, -3/5, -4/5])]
def pauli(v): return v[0] * X + v[1] * Y + v[2] * Z

def family(G):
    L = 4 * G - 20
    q = [(G - 1) / L] * 2 + [(G - 9) / L] * 2
    rho = [(I2 - pauli(s[i]) / (G - 1)) / 2 for i in range(2)] + [(I2 - 9 * pauli(s[i]) / (G - 9)) / 2 for i in range(2, 4)]
    M = [(I2 - pauli(s[i])) / 4 for i in range(4)]
    P = [(I2 + pauli(s[i])) / 2 for i in range(4)]
    return L, q, rho, M, P

def psd(A, tol=1e-12):
    return np.linalg.eigvalsh((A + A.conj().T) / 2).min() >= -tol

# A. exact rational check at G = 19 via sympy, plus numeric for other G
G = sp.Rational(19)
Xs, Ys, Zs = sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[0, -sp.I], [sp.I, 0]]), sp.Matrix([[1, 0], [0, -1]])
ss = [(sp.Rational(3, 5), 0, sp.Rational(4, 5)), (-sp.Rational(3, 5), 0, sp.Rational(4, 5)), (0, sp.Rational(3, 5), -sp.Rational(4, 5)), (0, -sp.Rational(3, 5), -sp.Rational(4, 5))]
Ps = lambda v: v[0] * Xs + v[1] * Ys + v[2] * Zs
L = 4 * G - 20
q = [(G - 1) / L] * 2 + [(G - 9) / L] * 2
rho = [(sp.eye(2) - Ps(ss[i]) / (G - 1)) / 2 for i in range(2)] + [(sp.eye(2) - 9 * Ps(ss[i]) / (G - 9)) / 2 for i in range(2, 4)]
M = [(sp.eye(2) - Ps(ss[i])) / 4 for i in range(4)]
P = [(sp.eye(2) + Ps(ss[i])) / 2 for i in range(4)]
K = G * sp.eye(2) / (2 * L)
r = [1 / L] * 2 + [9 / L] * 2
okA = sum(q) == 1 and sp.simplify(sum(M, sp.zeros(2)) - sp.eye(2)) == sp.zeros(2)
for i in range(4):
    okA &= sp.simplify(rho[i].trace()) == 1 and all(ev >= 0 for ev in rho[i].eigenvals())
    okA &= sp.simplify(K - q[i] * rho[i] - r[i] * P[i]) == sp.zeros(2)
    okA &= sp.simplify(P[i] * M[i]) == sp.zeros(2)
    okA &= sp.simplify(P[i] * P[i] - P[i]) == sp.zeros(2)
succ = sp.simplify(sum((q[i] * rho[i] * M[i]).trace() for i in range(4)))
okA &= succ == G / L and succ == sp.Rational(19, 56)
print('A exact G=19: states', [rho[i].tolist() for i in range(4)][:1], '... success', succ, 'ok' if okA else 'BAD')
if not okA: bad += 1

# B. Holevo-YKL test of M after random channels
rng = np.random.default_rng(20261012)
def ykl_pass(q, rho_out, M):
    Gam = sum(q[i] * rho_out[i] @ M[i] for i in range(4))
    if np.abs(Gam - Gam.conj().T).max() > 1e-10:
        return False
    return all(psd(Gam - q[i] * rho_out[i], 1e-10) for i in range(4))
def apply(K, A): return sum(k @ A @ k.conj().T for k in K)
def random_channel(scale):
    A = rng.normal(size=(8, 2)) + 1j * rng.normal(size=(8, 2))
    V = np.vstack([np.eye(2), np.zeros((6, 2))]) + scale * A        # isometry near identity Kraus
    Q, _ = np.linalg.qr(V)
    return [Q[2 * k:2 * k + 2, :] for k in range(4)]
for Gv in [19.0, 25.0, 60.0]:
    Lv, qv, rv, Mv, Pv = family(Gv)
    if not ykl_pass(qv, rv, Mv):
        bad += 1; print('identity fails YKL', Gv)
    passes = 0; trials = 0
    for scale in [1e-4, 1e-3, 1e-2, 1e-1, 1.0]:
        for _ in range(400):
            Kr = random_channel(scale); trials += 1
            if ykl_pass(qv, [apply(Kr, x) for x in rv], Mv):
                passes += 1
    named = []
    for p in [1e-3, 0.1, 0.5]:
        named.append([np.sqrt(1 - 3 * p / 4) * I2, np.sqrt(p / 4) * X, np.sqrt(p / 4) * Y, np.sqrt(p / 4) * Z])   # depolarizing
        named.append([np.sqrt(1 - p) * I2, np.sqrt(p) * Z])                                                     # dephasing
        named.append([np.array([[1, 0], [0, np.sqrt(1 - p)]]), np.array([[0, np.sqrt(p)], [0, 0]])])          # amplitude damping
        th = p; named.append([np.cos(th) * I2 - 1j * np.sin(th) * Z])                                            # unitary
    npass = sum(ykl_pass(qv, [apply(Kr, x) for x in rv], Mv) for Kr in named)
    print(f'B G={Gv}: identity passes; {passes}/{trials} random channels pass; {npass}/{len(named)} named channels pass')
    if passes or npass:
        bad += 1

# C. symbolic Choi quadratic forms (Lemma 3)
eps, tx, ty, tz = sp.symbols('epsilon t_x t_y t_z', real=True)
C = sp.Matrix([[2 - 9 * eps / 5 + 2 * tz, 2 * (tx - sp.I * ty), 0, 2 - 10 * eps],
               [2 * (tx + sp.I * ty), 9 * eps / 5 - 2 * tz, -8 * eps, 0],
               [0, -8 * eps, 9 * eps / 5, 0],
               [2 - 10 * eps, 0, 0, 2 - 9 * eps / 5]])
v = sp.Matrix([1, 0, 0, -1]); w = sp.Matrix([0, 1, 2, 0])
qv_ = sp.simplify((v.T * C * v)[0]); qw_ = sp.simplify((w.T * C * w)[0])
# Choi matrix of the affine Bloch map with T = diag-ish form and t, rebuilt independently
def choi_from_affine(T, t):
    Cm = sp.zeros(4, 4)
    basis = [sp.Matrix([[1, 0], [0, 0]]), sp.Matrix([[0, 1], [0, 0]]), sp.Matrix([[0, 0], [1, 0]]), sp.Matrix([[0, 0], [0, 1]])]
    def N(A):
        a0 = A.trace(); vec = [ (A * Xs).trace(), (A * Ys).trace(), (A * Zs).trace() ]
        out = sp.Matrix(T) * sp.Matrix(vec) + a0 * sp.Matrix(t)
        return (a0 * sp.eye(2) + out[0] * Xs + out[1] * Ys + out[2] * Zs) / 2
    for a in range(2):
        for b in range(2):
            E = sp.zeros(2); E[a, b] = 1
            blk = N(E)
            for c in range(2):
                for d in range(2):
                    Cm[2 * a + c, 2 * b + d] = blk[c, d]
    return 2 * Cm
T = sp.Matrix([[1 - 9 * eps, 0, tx], [0, 1 - eps, ty], [0, 0, 1 - 9 * eps / 5 + tz]])
Cc = sp.simplify(choi_from_affine(T, [tx, ty, tz]))
okC = sp.simplify(qv_ + qw_ + 33 * eps / 5) == 0 and sp.simplify(Cc - C) == sp.zeros(4)
print('C Choi matrix rebuilt from the affine form:', 'matches' if sp.simplify(Cc - C) == sp.zeros(4) else 'DIFFERS', '; vCv + wCw =', sp.simplify(qv_ + qw_))
if not okC:
    bad += 1
print('bad=', bad)
sys.exit(0 if bad == 0 else 1)
