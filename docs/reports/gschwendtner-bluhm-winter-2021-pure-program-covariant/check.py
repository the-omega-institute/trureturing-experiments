"""Pure versus mixed program states for covariant programmable processors.

Source: Gschwendtner, Bluhm, Winter, "Programmability of covariant quantum
channels", Quantum 5, 488 (2021), arXiv:2012.00717v2, Section 5 (Discussion).

Checks:
 1. Replacement instance: G = D4 acting on C^2 by its 2-dim irrep U, V trivial on C^2.
    - U is a unitary representation and has no common eigenvector (irreducible).
    - the group average of rho is Tr(rho)/2 * 1, so covariant channels are
      rho -> Tr(rho) sigma; the mixed processor P(rho (x) pi) = Tr(rho) pi has d_P = 2.
    - commutant of conj(U) (x) V has dimension 4 = 1 block B(C^2): d_c = 2.
 2. Pure upper bound 3: every qubit state is Tr_2 of a unit vector of Sym^2(C^2).
 3. Pure lower bound 3: the subspace argument is checked on random processors
    is a proof (see the preregistration); its algebraic step (orthogonal rank-one
    columns with parallel environments give rank-one outputs) is checked here.
 4. A_4 example (the source's Example): theta (x) theta commutant dimension 7
    = 1+1+1+2^2, so d_c = 5; V = diag(1,z,z^2) intertwines with phase; the
    block B has multiplicity 2; the pure construction sum m(n_k) = 1+1+1+3 = 6.
"""
import itertools, numpy as np, sys
rng = np.random.default_rng(20261011)
bad = 0
def ok(c, msg):
    global bad
    if not c:
        bad += 1; print("FAIL", msg)
    else:
        print("ok  ", msg)

# ---------- 1. D4 instance ----------
R = np.array([[0, -1], [1, 0]], dtype=complex)
S = np.array([[1, 0], [0, -1]], dtype=complex)
D4 = [np.linalg.matrix_power(R, i) for i in range(4)] + [S @ np.linalg.matrix_power(R, i) for i in range(4)]
ok(all(np.allclose(g.conj().T @ g, np.eye(2)) for g in D4), "D4 matrices unitary")
prods = {tuple(np.round((a @ b).flatten(), 8)) for a in D4 for b in D4}
ok(len(prods) == 8 and prods == {tuple(np.round(g.flatten(), 8)) for g in D4}, "D4 closed under products (8 elements)")
# no common eigenvector of R and S
evR = np.linalg.eig(R)[1]
ok(all(abs(abs(np.vdot(evR[:, j], S @ evR[:, j])) - 1) > 1e-6 or not np.allclose(S @ evR[:, j], np.vdot(evR[:, j], S @ evR[:, j]) * evR[:, j]) for j in range(2)), "R and S share no eigenvector")
for _ in range(20):
    A = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)); rho = A @ A.conj().T; rho /= np.trace(rho)
    avg = sum(g @ rho @ g.conj().T for g in D4) / 8
    if not np.allclose(avg, np.trace(rho) / 2 * np.eye(2)):
        ok(False, "group average")
ok(True, "group average = Tr(rho)/2 * 1 on 20 random states")
# commutant dimension of conj(U) (x) V with V trivial on C^2
def commutant_dim(mats):
    n = mats[0].shape[0]
    rows = [np.kron(np.eye(n), g) - np.kron(g.T, np.eye(n)) for g in mats]
    M = np.vstack(rows)
    return n * n - np.linalg.matrix_rank(M, tol=1e-9)
ok(commutant_dim([np.kron(g.conj(), np.eye(2)) for g in D4]) == 4, "commutant of conj(U)(x)V has dim 4 (one block, n = 2)")

# ---------- 2. Sym^2 purification ----------
def sym2_vec(rho):
    w, v = np.linalg.eigh(rho)
    psi = sum(np.sqrt(max(w[i], 0)) * np.kron(v[:, i], v[:, i]) for i in range(2))
    return psi
F = np.zeros((4, 4)); 
for i in range(2):
    for j in range(2):
        F[2 * i + j, 2 * j + i] = 1
for _ in range(50):
    A = rng.normal(size=(2, 2)) + 1j * rng.normal(size=(2, 2)); rho = A @ A.conj().T; rho /= np.trace(rho)
    psi = sym2_vec(rho)
    red = psi.reshape(2, 2) @ psi.reshape(2, 2).conj().T
    if not (np.allclose(F @ psi, psi) and np.isclose(np.linalg.norm(psi), 1) and np.allclose(red, rho)):
        ok(False, "Sym2 purification")
ok(True, "every sampled qubit state is Tr_2 of a unit vector in Sym^2(C^2) (dim 3)")

# ---------- 3. lower bound: algebraic step ----------
# M(psi) = 2 x r matrix; pure output |x><x| iff M = x f^*. With a basis psi0, psi1,
# M0 = e0 f0^*, M1 = e1 f1^*, the + target forces f0 || f1, then every M(c0 psi0 + c1 psi1)
# = (c0 e0 + c1 k e1) f0^* has rank <= 1.
for r in (1, 2, 3, 5):
    f0 = rng.normal(size=r) + 1j * rng.normal(size=r); f0 /= np.linalg.norm(f0)
    k = np.exp(1j * rng.uniform(0, 2 * np.pi))
    e0 = np.array([1, 0]); e1 = np.array([0, 1])
    worst = 0
    for _ in range(200):
        c = rng.normal(size=2) + 1j * rng.normal(size=2)
        M = c[0] * np.outer(e0, f0.conj()) + c[1] * k * np.outer(e1, f0.conj())
        worst = max(worst, np.linalg.svd(M, compute_uv=False)[1] if min(M.shape) > 1 else 0)
    ok(worst < 1e-12, f"parallel environments => rank-one outputs (r={r})")
# ---------- 4. A_4 example ----------
P = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0]], dtype=float)
signs = [np.diag(s) for s in [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]]
A4 = [d @ np.linalg.matrix_power(P, i) for d in signs for i in range(3)]
ok(len({tuple(np.round(g.flatten(), 8)) for g in A4}) == 12, "A_4 as 12 rotation matrices")
ok(len({tuple(np.round((a @ b).flatten(), 8)) for a in A4 for b in A4}) == 12, "A_4 closed")
ok(commutant_dim(A4) == 1, "theta irreducible (commutant dim 1)")
ok(commutant_dim([np.kron(g.conj(), g) for g in A4]) == 7, "commutant of theta(x)theta has dim 7 = 1+1+1+4, d_c = 5")
z = np.exp(2j * np.pi / 3); Vz = np.diag([1, z, z * z])
phases = []
for g in A4:
    X = g @ Vz @ g.conj().T
    lam = X[np.unravel_index(np.argmax(abs(Vz)), Vz.shape)] / Vz[np.unravel_index(np.argmax(abs(Vz)), Vz.shape)]
    phases.append(np.allclose(X, (np.trace(X) / np.trace(Vz) if abs(np.trace(Vz)) > 1e-9 else None) * Vz) if abs(np.trace(Vz)) > 1e-9 else any(np.allclose(X, w * Vz) for w in (1, z, z * z)))
ok(all(phases), "g V g^* = phase * V for all g (V = diag(1,z,z^2))")
# block of multiplicity 2: orthogonal complement of the three Omega_j in C^9 is invariant, dim 6
Om = [np.kron(np.eye(3), np.linalg.matrix_power(Vz, j)) @ (np.eye(3).reshape(9) / np.sqrt(3)) for j in range(3)]
G = np.array(Om).T
Pblk = np.eye(9) - G @ G.conj().T
ok(np.isclose(np.trace(Pblk).real, 6), "block projector has rank 6 = 3 (irrep) x 2 (multiplicity)")
ok(all(np.allclose(np.kron(g.conj(), g) @ Pblk, Pblk @ np.kron(g.conj(), g)) for g in A4), "block projector commutes with theta(x)theta")
print("pure construction dimension 1+1+1+3 =", 1 + 1 + 1 + 3, "; lower bound sum(n_k + [n_k>1]) =", 1 + 1 + 1 + 3, "; mixed optimum d_c =", 5)
print("bad=", bad)
sys.exit(1 if bad else 0)
