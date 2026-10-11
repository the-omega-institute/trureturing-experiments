#!/usr/bin/env python3
"""Chiribella, Davidson, Paulsen, Rahaman, "Positive maps and entanglement in real Hilbert spaces", arXiv:2207.02510v2;
Ann. Henri Poincare 24 (2023), Section 9.2: PPT^2(C) => IPT^2(C) => IPT^2(R) (Proposition 9.9), and "We do not know if
these conjectures are equivalent."

The flag lift. e_pm = (1, +-i)/sqrt2, P_pm = e_pm e_pm^*. For X on C^n (x) C^2, B_pm(X) = (I (x) e_pm^*) X (I (x) e_pm),
E(X) = B_+(X) + B_-(X)^T, D(A) = A (x) P_+ + A^T (x) P_-, Gamma_Phi = D Phi E / 2.

Checks (random instances, double precision, tolerance 1e-9):
  1. E D = 2 id, E tau = E, tau D = D, Tr(M E(X)) = Tr(D(M) X).
  2. For random PPT maps Phi (Choi and partial-transposed Choi positive): Gamma_Phi is completely positive (Choi PSD),
     input-transpose invariant (Gamma tau = Gamma, i.e. IPT), output symmetric, and maps real matrices to real matrices.
  3. Composition: Gamma_Psi Gamma_Phi = Gamma_(Psi Phi).
  4. Recovery: B_+ Gamma_Phi^r (. (x) P_+) = Phi^r / 2 for r = 1, 2, 3.
  5. Measure-and-prepare transfer: for an entanglement-breaking Phi(A) = sum Tr(M_j A) N_j, Gamma_Phi(X) =
     1/2 sum Tr(D(M_j) X) D(N_j) with D(M_j), D(N_j) real positive semidefinite.
  6. The source's idempotent PPT example before Conjecture IPT2 (n = 2, C = 1/2 [[I, g], [-g, I]]) is PPT and not IPT,
     and its flag lift is a real IPT map."""
import sys, itertools
import numpy as np

rng = np.random.default_rng(20261011)
bad = 0
TOL = 1e-9
def report(label, ok):
    global bad
    print(label, ok)
    bad += 0 if ok else 1

ep = np.array([1, 1j]) / np.sqrt(2); em = np.array([1, -1j]) / np.sqrt(2)
Pp = np.outer(ep, ep.conj()); Pm = np.outer(em, em.conj())

def Bc(X, n, e):            # (I_n (x) e^*) X (I_n (x) e), X on C^n (x) C^2 (flag index last)
    V = np.kron(np.eye(n), e.reshape(2, 1))          # 2n x n
    return V.conj().T @ X @ V
def E(X, n): return Bc(X, n, ep) + Bc(X, n, em).T
def D(A): return np.kron(A, Pp) + np.kron(A.T, Pm)

def choi(f, n, m):          # C = sum_ij E_ij (x) f(E_ij)
    C = np.zeros((n * m, n * m), dtype=complex)
    for i, j in itertools.product(range(n), repeat=2):
        Eij = np.zeros((n, n), dtype=complex); Eij[i, j] = 1
        C += np.kron(Eij, f(Eij))
    return C
def ptr_first(C, n, m):    # partial transpose on the first factor
    return C.reshape(n, m, n, m).transpose(2, 1, 0, 3).reshape(n * m, n * m)
psd = lambda C: np.allclose(C, C.conj().T, atol=TOL) and np.linalg.eigvalsh((C + C.conj().T) / 2).min() > -1e-8

def kraus_map(K):
    return lambda A: sum(k @ A @ k.conj().T for k in K)
def random_ppt(n):
    # mixture of a random CP map and the transpose-composed version of another: Phi = Phi1 + Phi2 o tau is PPT
    K1 = [rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)) for _ in range(2)]
    K2 = [rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n)) for _ in range(2)]
    f1, f2 = kraus_map(K1), kraus_map(K2)
    # Phi = f1 o tau o f1' style is not generally PPT; use Choi-level construction: C = S + S^Gamma with S separable-free PSD
    # Simpler PPT family: Phi(A) = f1(A) + f1(A^T) composed: Choi(Phi) = C1 + C1^Gamma, both PSD iff C1^Gamma PSD.
    # Use maps of the form A -> sum_j (Tr(M_j A) N_j) + t * f1(A) with small t; check PPT explicitly and resample.
    for _ in range(200):
        Ms = [(lambda G: G @ G.conj().T)(rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))) for _ in range(3)]
        Ns = [(lambda G: G @ G.conj().T)(rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))) for _ in range(3)]
        t = 0.05
        f = lambda A, Ms=Ms, Ns=Ns: sum(np.trace(M @ A) * N for M, N in zip(Ms, Ns)) + t * f1(A)
        C = choi(f, n, n)
        if psd(C) and psd(ptr_first(C, n, n)):
            return f
    raise RuntimeError('no PPT sample')

n = 3
# 1. identities
ok = True
for _ in range(50):
    A = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    X = rng.normal(size=(2 * n, 2 * n)) + 1j * rng.normal(size=(2 * n, 2 * n))
    M = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
    ok &= np.allclose(E(D(A), n), 2 * A) and np.allclose(E(X.T, n), E(X, n)) and np.allclose(D(A).T, D(A))
    ok &= abs(np.trace(M @ E(X, n)) - np.trace(D(M) @ X)) < 1e-9
report('1. E D = 2 id, E tau = E, tau D = D, Tr(M E X) = Tr(D(M) X):', ok)

# 2-4. lift properties
ok2 = ok3 = ok4 = True
for _ in range(10):
    Phi = random_ppt(n); Psi = random_ppt(n)
    G = lambda X, Phi=Phi: D(Phi(E(X, n))) / 2
    CG = choi(G, 2 * n, 2 * n)
    ok2 &= psd(CG) and np.allclose(ptr_first(CG, 2 * n, 2 * n), CG, atol=1e-8)
    for _ in range(5):
        X = rng.normal(size=(2 * n, 2 * n))
        Y = G(X)
        ok2 &= np.allclose(Y.imag, 0, atol=1e-9) and np.allclose(Y, Y.T, atol=1e-9)
    GPsi = lambda X, Psi=Psi: D(Psi(E(X, n))) / 2
    for _ in range(5):
        X = rng.normal(size=(2 * n, 2 * n)) + 1j * rng.normal(size=(2 * n, 2 * n))
        ok3 &= np.allclose(GPsi(G(X)), D(Psi(Phi(E(X, n)))) / 2, atol=1e-8)
    for r in (1, 2, 3):
        for _ in range(3):
            A = rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))
            Y = np.kron(A, Pp)
            for _ in range(r): Y = G(Y)
            Z = A.copy()
            for _ in range(r): Z = Phi(Z)
            ok4 &= np.allclose(Bc(Y, n, ep), Z / 2, atol=1e-6 * (1 + abs(Z).max()))
report('2. Gamma_Phi is CP, IPT (partial-transpose invariant Choi) and real-preserving with symmetric output:', ok2)
report('3. Gamma_Psi Gamma_Phi = Gamma_(Psi Phi):', ok3)
report('4. B_+ Gamma_Phi^r ( . (x) P_+) = Phi^r / 2 for r = 1, 2, 3:', ok4)

# 5. measure-and-prepare transfer
ok = True
for _ in range(10):
    Ms = [(lambda G_: G_ @ G_.conj().T)(rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))) for _ in range(4)]
    Ns = [(lambda G_: G_ @ G_.conj().T)(rng.normal(size=(n, n)) + 1j * rng.normal(size=(n, n))) for _ in range(4)]
    Phi = lambda A: sum(np.trace(M @ A) * N for M, N in zip(Ms, Ns))
    for _ in range(5):
        X = rng.normal(size=(2 * n, 2 * n)) + 1j * rng.normal(size=(2 * n, 2 * n))
        lhs = D(Phi(E(X, n))) / 2
        rhs = sum(np.trace(D(M) @ X) * D(N) for M, N in zip(Ms, Ns)) / 2
        ok &= np.allclose(lhs, rhs, atol=1e-8)
    ok &= all(np.allclose(D(M).imag, 0, atol=1e-12) and psd(D(M)) for M in Ms + Ns)
report('5. Gamma of a measure-and-prepare map is real measure-and-prepare with D(M_j), D(N_j) real PSD:', ok)

# 6. the source's real PPT idempotent example (n = 2)
g = np.array([[0, 1], [-1, 0]], dtype=float)
C = 0.5 * np.block([[np.eye(2), g], [-g, np.eye(2)]])
def from_choi(C, n, m):
    Cr = C.reshape(n, m, n, m)
    return lambda A: np.einsum('ij,iajb->ab', A, Cr)
Phi = from_choi(C.astype(complex), 2, 2)
ok = psd(choi(Phi, 2, 2)) and psd(ptr_first(choi(Phi, 2, 2), 2, 2)) and not np.allclose(ptr_first(C, 2, 2), C)
G = lambda X: D(Phi(E(X, 2))) / 2
CG = choi(G, 4, 4)
ok &= psd(CG) and np.allclose(ptr_first(CG, 4, 4), CG, atol=1e-9) and np.allclose(CG.imag, 0, atol=1e-12)
report('6. the source example map is PPT, not IPT, and its flag lift is a real IPT map:', ok)
print('bad=', bad)
sys.exit(1 if bad else 0)
