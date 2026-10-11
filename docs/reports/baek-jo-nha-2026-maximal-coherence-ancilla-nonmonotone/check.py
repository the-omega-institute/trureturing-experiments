#!/usr/bin/env python3
"""Baek, Jo, Nha, "Maximal coherence of quantum measurement and the resource theory of sharpness", arXiv:2607.05847v1,
Section III.C: "whether C_D^max remains monotonic under general CP unital maps in this regime is left as an open
problem" (Table I: non-dichotomic POVMs, any fuzzifying operation, "Open").

Distance D_inf(A, B) = 1/2 sum_x ||A_x - B_x||_op (eq. (DpSinglePOVM), p = infinity; Section IV.A uses it).
C(A) = min over incoherent (computational-basis diagonal) POVMs B of D_inf(A, B); C^max(A) = max_U C(U^dag A U);
S(A) = min over trivial POVMs {q_x I} of D_inf(A, {q_x I}).

Witness: the paper's tetrahedral SIC qubit POVM A_i = (I + v_i . sigma)/4 and the CP unital map J(X) = X (x) I_2
(Heisenberg dual of discarding a qubit; Kraus operators I (x) <j|).

Checks:
  1. A is a POVM; J is CP unital (Kraus form, sum K_j^dag K_j = I) and J(A) is a POVM on C^4.
  2. S(A) = 1/2: ||A_i - q I||_op = |1/4 - q| + 1/4 exactly, so S = 1/2 at q = 1/4 (exact, Fractions).
  3. Upper bound C^max(A) <= sqrt(2/3)/2: for any orthonormal a, b in C^2, sum_k |<a|sigma_k|b>|^2 = 2, so
     sum_i |<a|A_i|b>|^2 = 1/6; the diagonal part of U^dag A U is an incoherent POVM at distance
     1/2 sum_i |<a|A_i|b>| <= 1/2 * 2 * sqrt(1/6) = sqrt(2/3)/2 < 1/2. Numerically on 20000 random unitaries:
     the identity holds and the dephasing distance never exceeds the bound; the bound is attained (u = x axis).
  4. Lower bound C^max(J(A)) >= 1/2: in the Bell basis, every incoherent POVM B has R(B_x) = (Tr B_x / 4) I with
     R(Y) = Tr_2(Y)/2, and for every unit psi, ||Y||_op >= 1/2 sum_j <psi (x) j|Y|psi (x) j> = <psi|R(Y)|psi>.
     Checked: R of each Bell projector is I/4 (exact); on 20000 random Bell-diagonal POVMs the distance is >= 1/2;
     the trivial POVM {I/4} attains 1/2.
  5. Hence C^max(J(A)) >= 1/2 > sqrt(2/3)/2 >= C^max(A): J increases C^max, so C^max_{D_inf} is not monotone.
  6. Stabilization identity for random qubit POVMs: max over the Bell basis lower bound equals S(A) (sampled)."""
import sys, itertools, math
from fractions import Fraction as F
import numpy as np

rng = np.random.default_rng(20261011)
bad = 0
def report(label, ok):
    global bad
    print(label, ok)
    bad += 0 if ok else 1

I2 = np.eye(2, dtype=complex)
sx = np.array([[0, 1], [1, 0]], dtype=complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.array([[1, 0], [0, -1]], dtype=complex)
sig = [sx, sy, sz]
V = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]]) / math.sqrt(3)
A = [(I2 + sum(V[i, k] * sig[k] for k in range(3))) / 4 for i in range(4)]
opn = lambda X: np.linalg.norm(X, 2)
D = lambda P, Q: 0.5 * sum(opn(p - q) for p, q in zip(P, Q))

# 1. POVMs and the CP unital map
psd = lambda X: np.allclose(X, X.conj().T) and min(np.linalg.eigvalsh(X)) > -1e-12
J = lambda X: np.kron(X, I2)
kraus = [np.kron(I2, np.array([[1, 0]], dtype=complex)), np.kron(I2, np.array([[0, 1]], dtype=complex))]   # 2 x 4
ok = all(psd(a) for a in A) and np.allclose(sum(A), I2)
ok &= np.allclose(sum(K.conj().T @ K for K in kraus), np.eye(4))
ok &= all(np.allclose(sum(K.conj().T @ X @ K for K in kraus), J(X)) for X in A)
ok &= all(psd(J(a)) for a in A) and np.allclose(sum(J(a) for a in A), np.eye(4))
report('1. A is a POVM; J(X) = X (x) I_2 = sum_j K_j^dag X K_j with sum K_j^dag K_j = I; J(A) is a POVM:', ok)

# 2. exact sharpness: eigenvalues of A_i - qI are 1/4 - q +- 1/4, so the op norm is |1/4 - q| + 1/4
ok = all(abs(opn(A[i] - q * I2) - (abs(0.25 - q) + 0.25)) < 1e-12 for i in range(4) for q in np.linspace(0, 1, 41))
S_exact = min(F(1, 2) * sum(abs(F(1, 4) - q) + F(1, 4) for q in qs)
              for qs in [(F(1, 4),) * 4, (F(1, 2), F(1, 2), 0, 0), (1, 0, 0, 0)])
ok &= S_exact == F(1, 2)
report('2. ||A_i - qI|| = |1/4 - q| + 1/4 and S(A) = 1/2 (attained at q = 1/4; the sum of |1/4 - q_i| is >= 0):', ok)

# 3. upper bound on C^max(A)
def haar(n):
    Z = (rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))) / math.sqrt(2)
    Q, R = np.linalg.qr(Z); return Q * (np.diag(R) / abs(np.diag(R)))
bound = math.sqrt(2 / 3) / 2
ok = True; best = 0
for _ in range(20000):
    U = haar(2); a, b = U[:, 0], U[:, 1]
    ok &= abs(sum(abs(a.conj() @ s @ b) ** 2 for s in sig) - 2) < 1e-10
    ok &= abs(sum(abs(a.conj() @ Ai @ b) ** 2 for Ai in A) - 1 / 6) < 1e-10
    R = [U.conj().T @ Ai @ U for Ai in A]
    B = [np.diag(np.diag(r)) for r in R]
    ok &= psd(B[0]) and all(psd(x) for x in B) and np.allclose(sum(B), I2)
    dist = D(R, B); best = max(best, dist)
    ok &= dist <= bound + 1e-12
Ux = np.array([[1, 1], [1, -1]], dtype=complex) / math.sqrt(2)          # reference axis along x: (v_i . x)^2 = 1/3
Rx = [Ux.conj().T @ Ai @ Ux for Ai in A]
ok &= abs(D(Rx, [np.diag(np.diag(r)) for r in Rx]) - bound) < 1e-12
print(f'   sampled max dephasing distance {best:.6f}, bound sqrt(2/3)/2 = {bound:.6f}')
report('3. C^max(A) <= sqrt(2/3)/2 = 0.408248... (identity, POVM property of the dephased POVM, attainment):', ok)

# 4. lower bound on C^max(J(A)) in the Bell basis
bell = [np.array(v, dtype=complex) / math.sqrt(2) for v in ([1, 0, 0, 1], [1, 0, 0, -1], [0, 1, 1, 0], [0, 1, -1, 0])]
UB = np.column_stack(bell)
def ptr2(Y):                                     # R(Y) = Tr_2(Y)/2 on C^2 (x) C^2
    return np.einsum('ajbj->ab', Y.reshape(2, 2, 2, 2)) / 2
ok = all(np.allclose(ptr2(np.outer(b, b.conj())), I2 / 4) for b in bell)
JA = [UB.conj().T @ J(Ai) @ UB for Ai in A]         # U^dag J(A) U; incoherent B diagonal <-> Bell-diagonal U B U^dag
low = 1e9
for _ in range(20000):
    W = rng.dirichlet(np.ones(4), size=4).T          # W[x, k] >= 0, columns sum to 1: diag POVM entries
    B = [np.diag(W[x]).astype(complex) for x in range(4)]
    d_ = D(JA, B); low = min(low, d_)
    Y = [UB @ B[x] @ UB.conj().T for x in range(4)]
    ok &= all(np.allclose(ptr2(Y[x]), np.trace(B[x]) / 4 * I2) for x in range(4))
    ok &= d_ >= 0.5 - 1e-12
triv = [np.eye(4, dtype=complex) / 4] * 4
ok &= abs(D(JA, triv) - 0.5) < 1e-12
print(f'   sampled min over incoherent B of D(U^dag J(A) U, B): {low:.6f}')
report('4. C^max(J(A)) >= 1/2 (Bell basis; R maps Bell projectors to I/4; trivial POVM attains 1/2):', ok)

# 5. strict increase
report('5. C^max(J(A)) >= 1/2 > sqrt(2/3)/2 >= C^max(A): J increases maximal coherence:', 0.5 > bound)

# 6. stabilization on random qubit POVMs: the Bell-basis lower bound reaches S(A)
ok = True
for _ in range(200):
    G = [haar(2) @ np.diag(rng.uniform(0.1, 1, 2)) @ haar(2).conj().T for _ in range(4)]
    G = [g @ g.conj().T for g in G]; T = sum(G); w, Q = np.linalg.eigh(T); Tm = Q @ np.diag(w ** -0.5) @ Q.conj().T
    P = [Tm @ g @ Tm for g in G]
    S = min(0.5 * sum(opn(p - q * I2) for p, q in zip(P, qs)) for qs in rng.dirichlet(np.ones(4), size=3000))
    JP = [UB.conj().T @ J(p) @ UB for p in P]
    m = min(D(JP, [np.diag(W[x]).astype(complex) for x in range(4)]) for W in [rng.dirichlet(np.ones(4), size=4).T for _ in range(300)])
    ok &= m >= S - 1e-3
report('6. random qubit POVMs: Bell-basis coherence of P (x) I_2 is at least S(P) (sampled):', ok)
print('bad=', bad)
sys.exit(1 if bad else 0)
