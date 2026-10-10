#!/usr/bin/env python3
"""Han–Liu, arXiv:2311.16349v1 (Rep. Math. Phys. 93 (2024) 87–102), Section 5, final Conjecture:
for pi = m_1 pi_1 + ... + m_d pi_d (pairwise inequivalent irreducibles) with I_k <= d < I_{k+1},
pr(Phi_pi) = max{beta(Phi_pi), k}, where beta(Phi_pi) = max m_i (their Theorem 3.x) and I_k is the minimal length
of a phase-retrievable frame for C^k.

Counterexample family: G = C_3 with its three characters, each of multiplicity 2s. Then d = 3, k = 1, beta = 2s,
but pr >= 3s. More precisely, for the three-block pinching with blocks m_1, m_2, m_3 (N = sum, r = max),
pr = max{r, floor(N/2)}.

(1) The C_3 twirl Phi(X) = (1/3) sum_g pi(g) X pi(g)^* equals the three-block pinching (exact, sympy, small m).
(2) Lower bound: the encoding V(u,v,w) = ((u,v),(u,w),(v,w))/sqrt2 of C^a+C^b+C^c (a = n-m3', b = n-m2',
    c = n-m1', sum = n = floor(N/2)) is an isometry and X -> Phi(V X V^*) is injective on all operators
    (exact rank = n^2 over C) for every block triple with N <= 12 and n > r.
(3) Upper bound mechanism: for random subspaces M of dimension q > max{r, floor(N/2)} the two steps of the proof
    produce explicit pairs x != y (as rank-one states) with Phi(xx*) = Phi(yy*): the u +- v pair when M meets a
    block, otherwise x and Dx from M cap D^{-1}M.
(4) I_2 >= 4: for random triples of vectors in C^2, the kernel of the three intensity functionals contains an
    indefinite Hermitian K = xx* - yy*, giving a collision; with 4 generic vectors no collision exists (rank 4).
(5) beta <= max m_i on the instance m = (2,2,2): any 3-dimensional subspace violates the Knill–Laflamme condition
    P P_i P = alpha_i P (checked on random subspaces); the general statement is the source's theorem."""
import sys, itertools, random
import numpy as np
import sympy as sp

bad = 0
def report(label, ok):
    global bad
    print(label, ok)
    bad += 0 if ok else 1

# (1) C_3 twirl = pinching (numerical, random complex X)
rng1 = np.random.default_rng(1)
for m in [(1, 1, 1), (2, 1, 1), (2, 2, 2), (3, 1, 2)]:
    N = sum(m); lab = np.concatenate([np.full(mi, i) for i, mi in enumerate(m)])
    ok1 = True
    for _ in range(5):
        X = rng1.normal(size=(N, N)) + 1j * rng1.normal(size=(N, N))
        tw = sum(np.diag(np.exp(2j * np.pi * g * lab / 3)) @ X @ np.diag(np.exp(-2j * np.pi * g * lab / 3)) for g in range(3)) / 3
        pin = np.where(lab[:, None] == lab[None, :], X, 0)
        ok1 &= np.allclose(tw, pin)
    report(f'(1) m={m}: C_3 twirl equals the three-block pinching:', ok1)

# (2) lower bound encoding
def encoding(m):
    N = sum(m); n = N // 2; mm = list(m)
    if N % 2 == 1:
        i = max(range(3), key=lambda t: mm[t]); mm[i] -= 1
    a, b, c = n - mm[2], n - mm[1], n - mm[0]
    assert min(a, b, c) >= 0 and a + b + c == n
    # logical coordinates: u (a), v (b), w (c); block 1 gets (u, v), block 2 gets (u, w), block 3 gets (v, w)
    offs = [0, m[0], m[0] + m[1]]
    V = np.zeros((N, n))
    U = list(range(0, a)); Vv = list(range(a, a + b)); W = list(range(a + b, n))
    def place(block, parts):
        row = offs[block]
        for part in parts:
            for col in part:
                V[row, col] = 1.0; row += 1
    place(0, [U, Vv]); place(1, [U, W]); place(2, [Vv, W])
    return V / np.sqrt(2), n

def pinch(X, m):
    out = np.zeros_like(X); o = 0
    for mi in m:
        out[o:o + mi, o:o + mi] = X[o:o + mi, o:o + mi]; o += mi
    return out

checked = 0; ok2 = True
for m in itertools.product(range(1, 9), repeat=3):
    N = sum(m)
    if N > 12: continue
    n, r = N // 2, max(m)
    if n <= r: continue
    V, nn = encoding(m)
    ok = np.allclose(V.T @ V, np.eye(nn))
    cols = []
    for p in range(nn):
        for q in range(nn):
            E = np.zeros((nn, nn)); E[p, q] = 1
            cols.append(pinch(V @ E @ V.T, m).ravel())
    A = np.array(cols).T
    ok &= np.linalg.matrix_rank(A) == nn * nn
    if not ok: ok2 = False; print('(2) failure at', m)
    checked += 1
report(f'(2) encoding is an isometry and injective on all operators for all {checked} block triples with N <= 12 and floor(N/2) > max m:', ok2)
for s in (1, 2, 3):
    m = (2 * s, 2 * s, 2 * s); V, nn = encoding(m)
    print(f'(2) m={m}: logical dimension {nn} = 3s, beta = {2 * s}')

# (3) upper bound mechanism on random subspaces
rng = np.random.default_rng(5)
def rank_one(x): return np.outer(x, x.conj())
def collision_found(m, q):
    N = sum(m)
    Mb = rng.normal(size=(N, q)) + 1j * rng.normal(size=(N, q))
    Q, _ = np.linalg.qr(Mb)
    offs = [0, m[0], m[0] + m[1]]
    # step 1: M meets a block? (generic M does not; force one case by inserting a block vector half the time)
    lam = np.exp(2j * np.pi * np.array([0, 1, 2]) / 3)
    D = np.diag(np.concatenate([np.full(mi, lam[i]) for i, mi in enumerate(m)]))
    # M cap D^{-1} M: solve Q a = D^{-1} Q b
    A = np.hstack([Q, -np.linalg.inv(D) @ Q])
    _, s, vh = np.linalg.svd(A)
    null = vh[-1].conj()
    x = Q @ null[:q]
    if np.linalg.norm(x) < 1e-9: return False
    y = D @ x
    same_pinch = np.allclose(pinch(rank_one(x), m), pinch(rank_one(y), m))
    different = not np.allclose(rank_one(x), rank_one(y))
    in_M = np.allclose(Q @ (Q.conj().T @ y), y)
    return same_pinch and different and in_M
ok = True; tried = 0
for m in [(2, 2, 2), (3, 2, 2), (3, 3, 2), (4, 3, 3)]:
    N = sum(m); q = max(max(m), N // 2) + 1
    if 2 * q - N <= 0: continue
    for _ in range(20):
        tried += 1
        ok &= collision_found(m, q)
report(f'(3) {tried} random subspaces of dimension max(r, floor(N/2)) + 1: x, Dx collide under the pinching:', ok)

# (4) I_2 >= 4
def herm_basis2():
    return [np.array([[1, 0], [0, 0]]), np.array([[0, 0], [0, 1]]), np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]])]
ok3 = True
for _ in range(200):
    f = rng.normal(size=(3, 2)) + 1j * rng.normal(size=(3, 2))
    B = herm_basis2()
    A = np.array([[np.real(np.trace(Bk @ np.outer(fj, fj.conj()))) for Bk in B] for fj in f])
    _, s, vh = np.linalg.svd(A)
    K = sum(c * Bk for c, Bk in zip(vh[-1], B))
    ev = np.linalg.eigvalsh(K)
    ok3 &= ev[0] < -1e-9 and ev[1] > 1e-9
report('(4) three vectors in C^2: the intensity kernel contains an indefinite K = xx*-yy* (200 random triples):', ok3)
f4 = rng.normal(size=(4, 2)) + 1j * rng.normal(size=(4, 2))
A4 = np.array([[np.real(np.trace(Bk @ np.outer(fj, fj.conj()))) for Bk in herm_basis2()] for fj in f4])
report('(4) four generic vectors in C^2 give a rank-4 intensity map (no collision):', np.linalg.matrix_rank(A4) == 4)

# (5) beta <= 2 on m = (2,2,2): random 3-dim subspaces fail Knill-Laflamme
m = (2, 2, 2); N = 6
Ps = []
o = 0
for mi in m:
    P = np.zeros((N, N)); P[o:o + mi, o:o + mi] = np.eye(mi); Ps.append(P); o += mi
ok5 = True
for _ in range(200):
    Mb = rng.normal(size=(N, 3)) + 1j * rng.normal(size=(N, 3)); Q, _ = np.linalg.qr(Mb)
    fails = False
    for P in Ps:
        C = Q.conj().T @ P @ Q
        alpha = np.trace(C).real / 3
        if not np.allclose(C, alpha * np.eye(3), atol=1e-8): fails = True
    ok5 &= fails
report('(5) m=(2,2,2): 200 random 3-dimensional subspaces violate Knill–Laflamme:', ok5)
print('bad=', bad)
sys.exit(1 if bad else 0)
