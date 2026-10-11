#!/usr/bin/env python3
"""Stonner-Kysela-Weir-Novotny-Alber-Jex, Phys. Lett. A 384 (2020) 126179 (arXiv:1912.06014v1),
Sec. III.B: the RUO R2(X) = p H X H^* + (1-p) W X W^*, H = h(x)h, W = uv(x)uv, with
h e_k = exp(2^(k-d) pi i) e_k, u e_k = e_(k mod d + 1), v = A (+) I_(d-2), A in U(2) without zero entries.
Numerical check (not a proof) that ker(R2 - 1) = span{P_sym, P_asym} for even d, that every other
eigenvalue of R2 has modulus < 1 (so R2^n converges), including the case a f + b c = 0 of A.
Exit 0 and final line `bad= 0` when every case agrees."""
import numpy as np, sys
rng = np.random.default_rng(20261012)

def build(d, A, p):
    h = np.diag([np.exp(2.0 ** (k - d) * np.pi * 1j) for k in range(1, d + 1)])
    u = np.zeros((d, d), complex)
    for k in range(1, d + 1):
        u[(k % d), k - 1] = 1  # u e_k = e_{(k mod d)+1}, zero-based rows
    v = np.eye(d, dtype=complex); v[:2, :2] = A
    w = u @ v
    H, W = np.kron(h, h), np.kron(w, w)
    n = d * d
    # superoperator on vec(X) (column stacking): vec(M X M^*) = (conj(M) kron M) vec(X)
    S = p * np.kron(H.conj(), H) + (1 - p) * np.kron(W.conj(), W)
    F = np.zeros((n, n))
    for i in range(d):
        for j in range(d):
            F[i * d + j, j * d + i] = 1
    Ps, Pa = (np.eye(n) + F) / 2, (np.eye(n) - F) / 2
    return S, Ps, Pa

def random_A(balanced=False):
    phi = rng.uniform(0, 2 * np.pi)
    if balanced:
        r = 1 / np.sqrt(2)
        al, be = r * np.exp(1j * rng.uniform(0, 2 * np.pi)), r * np.exp(1j * rng.uniform(0, 2 * np.pi))
    else:
        t = rng.uniform(0.05, np.pi / 2 - 0.05)
        al, be = np.cos(t) * np.exp(1j * rng.uniform(0, 2 * np.pi)), np.sin(t) * np.exp(1j * rng.uniform(0, 2 * np.pi))
    return np.exp(1j * phi) * np.array([[al, be], [-np.conj(be), np.conj(al)]])

bad = 0; cases = 0
for d in [2, 3, 4, 5, 6, 8]:
    trials = 6 if d <= 6 else 2
    for t in range(trials):
        A = random_A(balanced=(t % 2 == 1))
        a, b, c, f = A[0, 0], A[0, 1], A[1, 0], A[1, 1]
        p = rng.uniform(0.05, 0.95)
        S, Ps, Pa = build(d, A, p)
        n = d * d
        vs, va = Ps.reshape(-1, order='F'), Pa.reshape(-1, order='F')
        ok_fix = np.allclose(S @ vs, vs, atol=1e-9) and np.allclose(S @ va, va, atol=1e-9)
        ev = np.linalg.eigvals(S)
        near1 = np.sum(np.abs(ev - 1) < 1e-7)
        peripheral = np.sum(np.abs(np.abs(ev) - 1) < 1e-7)
        gap = 1 - max(abs(x) for x in ev if abs(x - 1) >= 1e-7)
        sv = np.linalg.svd(S - np.eye(n * n), compute_uv=False)
        kernel_dim = int(np.sum(sv < 1e-8))
        cases += 1
        flag = ok_fix and near1 == 2 and peripheral == 2 and kernel_dim == 2 and gap > 1e-6
        if not flag:
            bad += 1
        print(f"d={d} balanced={t%2==1} |af+bc|={abs(a*f+b*c):.2e} p={p:.3f} kernel_dim={kernel_dim} eig1={near1} peripheral={peripheral} spectral_gap={gap:.3e} {'ok' if flag else 'BAD'}")
print('cases', cases)
print('bad=', bad)
sys.exit(0 if bad == 0 else 1)
