"""Exhaustive check for Jack, arXiv:2411.16284 (Semigroup Forum 112 (2026)), Problem 9.4:
"Is there a finite semigroup for which ~o is a congruence and S/~o is not cancelative?"
Otto conjugacy: a ~o b iff exists g, h in S with a g = g b and b h = h a.
For every labelled semigroup of order n <= 4 (all associative tables, backtracking with associativity pruning):
when ~o is a congruence, the quotient S/~o is an abelian group (hence cancellative).  Prints ALL_OK."""
import itertools
def semigroups(n):
    T = [[None]*n for _ in range(n)]
    cells = [(i, j) for i in range(n) for j in range(n)]
    def assoc_ok():
        for a in range(n):
            for b in range(n):
                ab = T[a][b]
                if ab is None: continue
                for c in range(n):
                    bc = T[b][c]
                    if bc is None: continue
                    l = T[ab][c]; r = T[a][bc]
                    if l is not None and r is not None and l != r: return False
        return True
    def bt(k):
        if k == len(cells):
            yield [row[:] for row in T]; return
        i, j = cells[k]
        for v in range(n):
            T[i][j] = v
            if assoc_ok(): yield from bt(k+1)
        T[i][j] = None
    yield from bt(0)
def otto(T, n):
    R = [[False]*n for _ in range(n)]
    for a in range(n):
        for b in range(n):
            R[a][b] = any(T[a][g] == T[g][b] for g in range(n)) and any(T[b][h] == T[h][a] for h in range(n))
    return R
ok = True; counts = {}
for n in range(1, 5):
    tot = cong = 0
    for T in semigroups(n):
        tot += 1
        R = otto(T, n)
        congruence = all(R[T[a][c]][T[b][d]] for a in range(n) for b in range(n) if R[a][b] for c in range(n) for d in range(n) if R[c][d])
        if not congruence: continue
        cong += 1
        # equivalence check and quotient
        assert all(R[a][a] for a in range(n))
        cls = {}
        for a in range(n):
            key = frozenset(b for b in range(n) if R[a][b]); cls.setdefault(key, a)
        reps = list(cls.values()); m = len(reps)
        idx = {key: k for k, key in enumerate(cls)}
        cl = lambda a: idx[frozenset(b for b in range(n) if R[a][b])]
        Q = [[cl(T[reps[x]][reps[y]]) for y in range(m)] for x in range(m)]
        # abelian group: commutative, identity, inverses
        comm = all(Q[x][y] == Q[y][x] for x in range(m) for y in range(m))
        ids = [e for e in range(m) if all(Q[e][x] == x == Q[x][e] for x in range(m))]
        grp = bool(ids) and all(any(Q[x][y] == ids[0] for y in range(m)) for x in range(m))
        canc = all(len({Q[x][y] for y in range(m)}) == m and len({Q[y][x] for y in range(m)}) == m for x in range(m))
        if not (comm and grp and canc): ok = False; print("COUNTEREXAMPLE", n, T)
    counts[n] = (tot, cong)
    print(f"order {n}: {tot} labelled semigroups, {cong} with ~o a congruence; every such quotient is an abelian group")
print("ALL_OK" if ok else "FAIL")
