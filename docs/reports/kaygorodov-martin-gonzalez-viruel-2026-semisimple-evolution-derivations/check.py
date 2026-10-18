"""Kaygorodov–Martín González–Viruel arXiv:2610.22248, Conjecture 5: numerical check of the proof's formula.

For an evolution algebra with natural basis e_1..e_n and squares e_i^2 = w_{c(i)} (classes c of
proportional squares, rescaled to equal squares), Der(A) is computed exactly as the null space of the
derivation equations (floating point, rank by SVD). When the Killing form of Der(A) is nondegenerate
(semisimple, char 0), the proof predicts
   Der(A) = (+)_a h_a,  h_a = {X in so(V_a) : X W_a = 0},  W_a = span{P_a w_b},
   h_a = so(m_a) with m_a = dim W_a^perp - dim(W_a cap W_a^perp)  (form: V_a orthonormal),
and every nonzero summand has m_a >= 3. Checked: dim Der = sum_a m_a(m_a-1)/2 and no m_a = 2.
Also checked: every sampled algebra with a zero square (Ann != 0) has non-semisimple Der.
"""
import itertools, math, sys
import numpy as np

rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 7)
TOL = 1e-8


def null_space(M):
    if M.size == 0:
        return np.eye(M.shape[1])
    u, s, vt = np.linalg.svd(M)
    r = int((s > TOL * max(1.0, s.max() if s.size else 1)).sum())
    return vt[r:].T


def derivations(n, cls, W):
    """cls[i] = class of e_i (or -1 if e_i^2 = 0); W[a] = square vector of class a (length n)."""
    sq = [W[cls[i]] if cls[i] >= 0 else np.zeros(n) for i in range(n)]
    rows = []
    # unknown D as n*n vector, D[j,i] = coefficient of e_j in D(e_i); index j*n+i
    def var(j, i): return j * n + i
    for i in range(n):
        for j in range(n):
            if i == j:
                # D(e_i^2) = 2 d_ii e_i^2 : sum_k (e_i^2)_k D(e_k) - 2 d_ii e_i^2 = 0
                for t in range(n):
                    row = np.zeros(n * n)
                    for kk in range(n):
                        row[var(t, kk)] += sq[i][kk]
                    row[var(i, i)] -= 2 * sq[i][t]
                    rows.append(row)
            elif i < j:
                # D(e_i) e_j + e_i D(e_j) = d_ji e_j^2 + d_ij e_i^2 = 0
                for t in range(n):
                    row = np.zeros(n * n)
                    row[var(j, i)] += sq[j][t]
                    row[var(i, j)] += sq[i][t]
                    rows.append(row)
    N = null_space(np.array(rows))
    return [N[:, c].reshape(n, n) for c in range(N.shape[1])]


def killing_nondegenerate(basis):
    d = len(basis)
    if d == 0:
        return True
    B = np.array([b.ravel() for b in basis]).T          # n^2 x d
    pinv = np.linalg.pinv(B)
    ad = []
    for x in basis:
        cols = []
        for y in basis:
            c = x @ y - y @ x
            coef = pinv @ c.ravel()
            if np.linalg.norm(B @ coef - c.ravel()) > 1e-6:
                raise RuntimeError("Der not closed under bracket")
            cols.append(coef)
        ad.append(np.array(cols).T)
    K = np.array([[np.trace(a @ b) for b in ad] for a in ad])
    s = np.linalg.svd(K, compute_uv=False)
    return s.min() > 1e-7 * max(1.0, s.max())


def predicted(n, cls, W):
    classes = sorted(set(c for c in cls if c >= 0))
    ms = []
    for a in classes:
        idx = [i for i in range(n) if cls[i] == a]
        P = np.array([[W[b][i] for i in idx] for b in classes])   # rows: P_a w_b in V_a coords
        Wa = P if P.size else np.zeros((0, len(idx)))
        # W_a^perp in V_a w.r.t. the standard form: null space of Wa
        perp = null_space(Wa) if Wa.shape[0] else np.eye(len(idx))
        dperp = perp.shape[1]
        # R = W_a cap W_a^perp: vectors in perp lying in rowspace(Wa)
        if Wa.shape[0] and dperp:
            rk_w = np.linalg.matrix_rank(Wa, tol=TOL)
            rk_sum = np.linalg.matrix_rank(np.vstack([Wa, perp.T]), tol=TOL)
            dR = rk_w + dperp - rk_sum
        else:
            dR = 0
        ms.append(dperp - dR)
    return ms


def sample():
    n = int(rng.integers(3, 9))
    ncls = int(rng.integers(1, n + 1))
    cls = [int(rng.integers(0, ncls)) for _ in range(n)]
    if rng.random() < 0.25:
        cls[int(rng.integers(0, n))] = -1                  # a zero square (Ann != 0)
    used = sorted(set(c for c in cls if c >= 0))
    remap = {c: k for k, c in enumerate(used)}
    cls = [remap[c] if c >= 0 else -1 for c in cls]
    W = []
    for a in range(len(used)):
        w = np.zeros(n)
        support = rng.choice(n, size=int(rng.integers(1, min(3, n) + 1)), replace=False)
        w[support] = rng.integers(-2, 3, size=len(support))
        if not w.any():
            w[support[0]] = 1
        W.append(w)
    # classes must have pairwise non-proportional squares
    for a, b in itertools.combinations(range(len(W)), 2):
        if np.linalg.matrix_rank(np.vstack([W[a], W[b]]), tol=TOL) < 2:
            return None
    return n, cls, W


def main():
    stats = dict(samples=0, semisimple=0, nonzero=0, ok=0, zero_sq_semisimple=0)
    hist = {}
    fails = []
    trials = int(sys.argv[2]) if len(sys.argv) > 2 else 4000
    for _ in range(trials):
        s = sample()
        if s is None:
            continue
        n, cls, W = s
        stats["samples"] += 1
        D = derivations(n, cls, W)
        if not killing_nondegenerate(D):
            continue
        stats["semisimple"] += 1
        if any(c < 0 for c in cls):
            stats["zero_sq_semisimple"] += 1 if len(D) > 0 else 0
            if len(D) > 0:
                fails.append(("annihilator", n, cls))
            continue
        ms = predicted(n, cls, W)
        if len(D) > 0:
            stats["nonzero"] += 1
        good = sum(m * (m - 1) // 2 for m in ms) == len(D) and all(m != 2 for m in ms)
        stats["ok"] += good
        if len(D) > 0:
            key = tuple(sorted((m for m in ms if m >= 3), reverse=True))
            hist[key] = hist.get(key, 0) + 1
        if not good:
            fails.append(("formula", n, cls, [w.tolist() for w in W], ms, len(D)))
    print(stats)
    print("nonzero semisimple Der by summand sizes m_a:", dict(sorted(hist.items(), key=lambda kv: -kv[1])))
    # Proposition 14: e_i^2 = u (i <= n), u^2 = u: all squares equal, one class; Der = so_n
    for m in range(1, 7):
        n = m + 1
        W0 = np.zeros(n); W0[m] = 1
        cls = [0] * (m + 1)
        D = derivations(n, cls, [W0])
        print(f"Proposition 14, n={m}: dim Der = {len(D)}, expected {m*(m-1)//2}, predicted m = {predicted(n, cls, [W0])}")
        if len(D) != m * (m - 1) // 2:
            fails.append(("prop14", m))
    print("RESULT", "PASS" if not fails else f"FAIL {fails[:3]}")
    sys.exit(0 if not fails else 1)


if __name__ == "__main__":
    main()
