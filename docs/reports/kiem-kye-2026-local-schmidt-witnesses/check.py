"""Numerical check for Kiem-Kye, arXiv:2608.14199v1, Section 5 conjecture.

BP_k: trace-one Hermitian W on C^m (x) C^n with <x, W x> >= 0 for Schmidt rank(x) <= k.
For rho = P_{E-perp}/dim E-perp on the boundary of BP_k, the conjecture asks that every
neighbourhood of rho meets BP_l \\ BP_{l+1} (k <= l < min(m, n)).

Construction tested (proposed proof): v in ker rho with SR(v) = r <= k, v_t = (v + t z)/|.|
with SR(v_t) = l + 1, A_t = (I - P_t) rho (I - P_t), Q_t = alpha_t I - P_t with
alpha_t = max_{SR(x) <= l, |x|=1} |<v_t, x>|^2 (= sum of the top-l squared Schmidt
coefficients of v_t), W_t = (A_t + t Q_t)/tr(A_t + t Q_t).
Prediction registered before the run: for each case, (a) a numerical maximization of
|<v_t, x>|^2 over SR <= l matches the top-l formula, (b) <v_t, W_t v_t> < 0, (c) the minimum
of <x, W_t x> over SR(x) <= l found by local search is >= -1e-9, (d) |W_t - rho| -> 0.

Usage: python3 check.py [seed [restarts steps]]   (numpy; overlap search default 20 restarts x 300 steps)
"""
import sys
import numpy as np


def schmidt(v, m, n):
    return np.linalg.svd(v.reshape(m, n), compute_uv=False)


def min_form_rank(W, m, n, l, rng, restarts=30, steps=400):
    """Local search for min <x, W x> over unit x of Schmidt rank <= l (x = A B^T, A m x l, B n x l)."""
    best = np.inf
    for _ in range(restarts):
        A = rng.standard_normal((m, l)) + 1j * rng.standard_normal((m, l))
        B = rng.standard_normal((n, l)) + 1j * rng.standard_normal((n, l))

        def val(A, B):
            x = (A @ B.T).reshape(-1)
            x = x / np.linalg.norm(x)
            return float(np.real(np.vdot(x, W @ x)))
        cur = val(A, B)
        s = 0.3
        for _ in range(steps):
            nA = A + s * (rng.standard_normal(A.shape) + 1j * rng.standard_normal(A.shape))
            nB = B + s * (rng.standard_normal(B.shape) + 1j * rng.standard_normal(B.shape))
            nv = val(nA, nB)
            if nv < cur:
                A, B, cur = nA, nB, nv
            else:
                s *= 0.99
        best = min(best, cur)
    return best


def max_overlap_rank(v, m, n, l, rng, restarts=None, steps=None):
    restarts = restarts or SEARCH[0]; steps = steps or SEARCH[1]
    best = 0.0
    for _ in range(restarts):
        A = rng.standard_normal((m, l)) + 1j * rng.standard_normal((m, l))
        B = rng.standard_normal((n, l)) + 1j * rng.standard_normal((n, l))

        def val(A, B):
            x = (A @ B.T).reshape(-1)
            x = x / np.linalg.norm(x)
            return abs(np.vdot(v, x)) ** 2
        cur = val(A, B)
        s = 0.3
        for _ in range(steps):
            nA = A + s * (rng.standard_normal(A.shape) + 1j * rng.standard_normal(A.shape))
            nB = B + s * (rng.standard_normal(B.shape) + 1j * rng.standard_normal(B.shape))
            nv = val(nA, nB)
            if nv > cur:
                A, B, cur = nA, nB, nv
            else:
                s *= 0.99
        best = max(best, cur)
    return best


def kron_vec(a, b):
    return np.kron(a, b)


def case(m, n, Ebasis, k, l, rng, label):
    N = m * n
    E = np.array(Ebasis).T
    Q, _ = np.linalg.qr(E)
    Pperp = np.eye(N) - Q @ Q.conj().T
    rho = Pperp / np.trace(Pperp).real
    # a kernel vector of Schmidt rank <= k (taken from E)
    v = Ebasis[0] / np.linalg.norm(Ebasis[0])
    r = int(np.sum(schmidt(v, m, n) > 1e-9))
    assert r <= k
    U, s, Vh = np.linalg.svd(v.reshape(m, n))
    q = l + 1 - r
    z = sum(np.kron(U[:, i], Vh[i, :]) for i in range(r, l + 1)) / np.sqrt(q)
    ok = True
    print(f"{label}: m={m} n={n} dim E={len(Ebasis)} SR(v)={r} k={k} l={l}")
    for t in (0.2, 0.05, 0.01):
        vt = (v + t * z) / np.sqrt(1 + t * t)
        sv = schmidt(vt, m, n)
        assert int(np.sum(sv > 1e-12)) == l + 1
        alpha = float(np.sum(np.sort(sv ** 2)[::-1][:l]))
        alpha_num = max_overlap_rank(vt, m, n, l, rng)
        P = np.outer(vt, vt.conj())
        R = np.eye(N) - P
        A = R @ rho @ R
        Qt = alpha * np.eye(N) - P
        M = A + t * Qt
        W = M / np.trace(M).real
        neg = float(np.real(np.vdot(vt, W @ vt)))
        mn = min_form_rank(W, m, n, l, rng)
        dist = float(np.max(np.abs(W - rho)))
        good = abs(alpha - alpha_num) < 1e-4 and neg < 0 and mn > -1e-9 and np.allclose(W, W.conj().T)
        ok &= good
        print(f"   t={t}: alpha={alpha:.6f} (search {alpha_num:.6f})  <v_t,W v_t>={neg:.3e}  "
              f"min_SR<=l <x,Wx>={mn:.3e}  max|W-rho|={dist:.3e}  {'ok' if good else 'FAIL'}")
    return ok


SEARCH = (20, 300)


def main():
    global SEARCH
    if len(sys.argv) > 2:
        SEARCH = (int(sys.argv[2]), int(sys.argv[3]))
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 20261016
    rng = np.random.default_rng(seed)
    e = np.eye(3)
    ok = True
    # non-generic: E spanned by one product vector (k = 1), 3 x 3
    prod = kron_vec(e[0], e[0])
    ok &= case(3, 3, [prod], 1, 1, rng, "product kernel")
    ok &= case(3, 3, [prod], 1, 2, rng, "product kernel")
    # E spanned by a product vector and an entangled vector, 3 x 3
    ent = (kron_vec(e[1], e[2]) + kron_vec(e[2], e[1])) / np.sqrt(2)
    ok &= case(3, 3, [prod, ent], 1, 2, rng, "mixed kernel")
    # rank-2 kernel vector only (k = 2), 3 x 4
    e4 = np.eye(4)
    v2 = (kron_vec(e[0], e4[0]) + 2 * kron_vec(e[1], e4[1])) / np.sqrt(5)
    ok &= case(3, 4, [v2], 2, 2, rng, "Schmidt-rank-2 kernel")
    print("ALL_OK" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
