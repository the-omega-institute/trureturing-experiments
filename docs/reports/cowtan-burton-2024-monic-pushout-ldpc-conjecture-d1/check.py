#!/usr/bin/env python3
"""Cowtan-Burton, arXiv:2301.13738v6, Conjecture D.1: the basis-preserving pushout of a monic span of LDPC Z-type
complexes is LDPC. Checks, on random based length-2 complexes over F2 and random basis-preserving monic spans C <- A -> D:
(1) the pushout differential dQ[y] = [d y] does not depend on the representative of a class (the chain-map condition);
(2) dQ1 dQ2 = 0; (3) column weight of dQ <= max(column weights of dC, dD) and row weight <= rowC + rowD, in each degree;
(4) the non-monic pushout of Appendix D (m pendant edges glued along m vertices to one vertex) has check weight m + 1.
Pure Python, exact F2 arithmetic. usage: check.py [trials]"""
import random, sys

def matmul(A, B, n, k, m):  # A: n x k, B: k x m, dict-of-sets rows -> F2
    return [[sum(A[i][t] & B[t][j] for t in range(k)) % 2 for j in range(m)] for i in range(n)]

def zeros(n, m): return [[0] * m for _ in range(n)]

def kernel(M, n, m):  # basis of {x in F2^m : M x = 0}
    rows = [r[:] for r in M]; piv = []; r = 0
    for c in range(m):
        p = next((i for i in range(r, n) if rows[i][c]), None)
        if p is None: continue
        rows[r], rows[p] = rows[p], rows[r]
        for i in range(n):
            if i != r and rows[i][c]: rows[i] = [(a ^ b) for a, b in zip(rows[i], rows[r])]
        piv.append(c); r += 1
    free = [c for c in range(m) if c not in piv]; basis = []
    for fc in free:
        x = [0] * m; x[fc] = 1
        for i, pc in enumerate(piv):
            if rows[i][fc]: x[pc] = 1
        basis.append(x)
    return basis

def rand_sparse_col(n, w, rng):
    c = [0] * n
    for i in rng.sample(range(n), min(w, n)): c[i] = 1
    return c

def rand_kernel_col(M, n, m, rng, w):
    K = kernel(M, n, m)
    if not K: return [0] * m
    x = [0] * m
    for v in rng.sample(K, min(w, len(K))): x = [a ^ b for a, b in zip(x, v)]
    return x

def make_A(rng):
    a = [rng.randint(0, 4) for _ in range(3)]  # dims a0 a1 a2
    d2 = [[0] * a[2] for _ in range(a[1])]
    for j in range(a[2]):
        col = rand_sparse_col(a[1], 2, rng)
        for i in range(a[1]): d2[i][j] = col[i]
    # d1 rows in the left kernel of d2: x^T d2 = 0  <=>  d2^T x = 0
    d2T = [[d2[i][j] for i in range(a[1])] for j in range(a[2])]
    d1 = [rand_kernel_col(d2T, a[2], a[1], rng, 2) for _ in range(a[0])]
    return a, d2, d1

def extend(a, d2A, d1A, rng):
    """A complex E containing A as a closed subcomplex, with extra basis elements, and the basis embedding of A,
    after a random relabelling of every basis."""
    x = [rng.randint(0, 4) for _ in range(3)]; e = [a[i] + x[i] for i in range(3)]
    d1 = zeros(e[0], e[1]); d2 = zeros(e[1], e[2])
    for i in range(a[0]):
        for j in range(a[1]): d1[i][j] = d1A[i][j]
    for j in range(a[1], e[1]):
        col = rand_sparse_col(e[0], 2, rng)
        for i in range(e[0]): d1[i][j] = col[i]
    for i in range(a[1]):
        for j in range(a[2]): d2[i][j] = d2A[i][j]
    for j in range(a[2], e[2]):
        col = rand_kernel_col(d1, e[0], e[1], rng, 2)
        for i in range(e[1]): d2[i][j] = col[i]
    perms = [rng.sample(range(e[i]), e[i]) for i in range(3)]  # new label of old index
    P1 = zeros(e[1], e[2]); Q1 = zeros(e[0], e[1])
    for i in range(e[1]):
        for j in range(e[2]): P1[perms[1][i]][perms[2][j]] = d2[i][j]
    for i in range(e[0]):
        for j in range(e[1]): Q1[perms[0][i]][perms[1][j]] = d1[i][j]
    emb = [[perms[k][i] for i in range(a[k])] for k in range(3)]
    return e, P1, Q1, emb

def row_w(M): return max([sum(r) for r in M], default=0)
def col_w(M, n, m): return max([sum(M[i][j] for i in range(n)) for j in range(m)], default=0)

def pushout(c, dC, d, dD, f, g):
    """Classes of (C ⊔ D)_k under f(a) ~ g(a); dQ[Y] = image of the column of a representative (both representatives
    are compared)."""
    cls = []
    for k in range(3):
        parent = {('C', i): ('C', i) for i in range(c[k])}; parent.update({('D', i): ('D', i) for i in range(d[k])})
        def find(u):
            while parent[u] != u: u = parent[u]
            return u
        for t in range(len(f[k])): parent[find(('D', g[k][t]))] = find(('C', f[k][t]))
        reps = {}
        for u in list(parent): reps.setdefault(find(u), []).append(u)
        cls.append(list(reps.values()))
    def col_of(u, deg):  # column of basis element u (degree deg) of dC or dD, as a set of (side, idx)
        side, i = u; M, n = (dC[deg], c[deg - 1]) if side == 'C' else (dD[deg], d[deg - 1])
        return {(side, r) for r in range(n) if M[r][i]}
    dQ = {}
    ok_welldef = True
    for deg in (1, 2):
        idx = {u: X for X, members in enumerate(cls[deg - 1]) for u in members}
        M = zeros(len(cls[deg - 1]), len(cls[deg]))
        for Y, members in enumerate(cls[deg]):
            imgs = []
            for u in members:
                v = [0] * len(cls[deg - 1])
                for w in col_of(u, deg): v[idx[w]] ^= 1
                imgs.append(v)
            if any(v != imgs[0] for v in imgs): ok_welldef = False
            for X in range(len(cls[deg - 1])): M[X][Y] = imgs[0][X]
        dQ[deg] = M
    return cls, dQ, ok_welldef

def trial(rng):
    a, d2A, d1A = make_A(rng)
    c, dC2, dC1, f = extend(a, d2A, d1A, rng)
    d, dD2, dD1, g = extend(a, d2A, d1A, rng)
    cls, dQ, wd = pushout(c, {2: dC2, 1: dC1}, d, {2: dD2, 1: dD1}, f, g)
    q = [len(cls[k]) for k in range(3)]
    zero = all(v == 0 for r in matmul(dQ[1], dQ[2], q[0], q[1], q[2]) for v in r)
    ok = wd and zero
    for deg, MC, MD, nc, nd in ((2, dC2, dD2, (c[1], c[2]), (d[1], d[2])), (1, dC1, dD1, (c[0], c[1]), (d[0], d[1]))):
        MQ = dQ[deg]; nq = (q[deg - 1], q[deg])
        ok &= col_w(MQ, *nq) <= max(col_w(MC, *nc), col_w(MD, *nd))
        ok &= row_w(MQ) <= row_w(MC) + row_w(MD)
    return ok, wd, zero

def appendix_d(m):
    # graph C: m edges, edge i joins vertices (2i, 2i+1); D the same; apex: m vertices, mapped to the vertices 2i in C and
    # ALL to vertex 0 in D (not injective), so the m edges of C and one edge of D meet one glued vertex: weight m + 1.
    # As complexes C1 -> C0 (edges -> vertices), with only degrees 1, 0 used.
    nC0, nC1 = 2 * m, m
    dC = [[1 if (v == 2 * e or v == 2 * e + 1) else 0 for e in range(nC1)] for v in range(nC0)]
    parent = {('C', v): ('C', v) for v in range(nC0)}; parent.update({('D', v): ('D', v) for v in range(nC0)})
    def find(u):
        while parent[u] != u: u = parent[u]
        return u
    for t in range(m): parent[find(('D', 0))] = find(('C', 2 * t))
    classes = {}
    for u in parent: classes.setdefault(find(u), []).append(u)
    deg_of = {X: 0 for X in classes}
    for side in ('C', 'D'):
        for e in range(nC1):
            for v in range(nC0):
                if dC[v][e]: deg_of[find((side, v))] += 1
    return max(deg_of.values()), row_w(dC)

if __name__ == "__main__":
    T = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
    rng = random.Random(20261020); bad = 0; nwd = 0; nz = 0
    for _ in range(T):
        ok, wd, z = trial(rng); bad += (not ok); nwd += wd; nz += z
    print(f"monic spans: {T} trials, well-defined {nwd}/{T}, dQ1 dQ2 = 0 in {nz}/{T}, bounds violated in {bad}")
    for m in (2, 3, 5, 10):
        wq, wc = appendix_d(m)
        print(f"Appendix D non-monic span, m = {m}: inputs have column weight 2 and vertex degree {wc}; glued vertex degree {wq} (= m + 1)")
    ok = bad == 0 and nwd == T and nz == T
    print("RESULT PASS" if ok else "RESULT FAIL"); print("EXIT=0" if ok else "EXIT=1"); sys.exit(0 if ok else 1)
