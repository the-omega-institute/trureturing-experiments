#!/usr/bin/env python3
"""Check of the proposed refutation of the ColorSym optimality conjecture of
Godin-Grisot-Togni (arXiv:2209.03436v1, Sec. 3.1; EJGTA 13(1) (2025) 57-74).

ColorSym (as written): Step 1 gives every vertex all colors of its singleton proper
intersection I_p({i}). Step i >= 2: if w_j >= b for every j, or i = n+1, stop; otherwise for
each set S with |S| = i, repeatedly remove a color of I_p(S) and give it to the j in S with
minimum w_j (smallest index on ties), until I_p(S) is empty.

Checks:
 1. The n=4, a=6, b=4 witness is an a-list assignment, is (L,b)-colorable (explicit coloring
    and brute-force max-flow), and ColorSym ends with vertex 2 holding 3 < b colors, for every
    order of same-size supports and every choice of removed color (counts only depend on the
    atom sizes here, since each support size has at most one nonempty atom).
 2. The family L_1=P+T, L_2=P+Q2, L_3=T+Q3, L_4=T+Q4 (a=3m, b=2m) for m=2..12: colorable,
    and ColorSym gives vertex 2 exactly m+floor(m/2) < 2m colors.
"""
import itertools, sys

def atoms(L):
    n = len(L); allc = set().union(*L)
    at = {}
    for c in sorted(allc):
        S = frozenset(i for i in range(n) if c in L[i])
        at.setdefault(S, []).append(c)
    return at

def colorsym(L, b, order_key=None, pick=min):
    n = len(L); at = atoms(L)
    w = [0]*n; out = [set() for _ in range(n)]
    for S, cs in at.items():
        if len(S) == 1:
            (i,) = tuple(S); out[i] |= set(cs); w[i] += len(cs)
    for size in range(2, n+1):
        sets = [S for S in at if len(S) == size]
        sets.sort(key=order_key or (lambda S: sorted(S)))
        for S in sets:
            pool = list(at[S])
            while pool:
                if all(x >= b for x in w):
                    return w, out
                c = pick(pool); pool.remove(c)
                j = min(S, key=lambda k: (w[k], k))
                out[j].add(c); w[j] += 1
    return w, out

def colorable(L, b):
    """Max-flow: source->vertex cap b, vertex->color if color in list, color->sink cap 1."""
    n = len(L); colors = sorted(set().union(*L))
    # simple augmenting paths on bipartite multi-matching
    cap_v = [b]*n; match = {}  # color -> vertex
    def try_assign(v, seen):
        for c in sorted(L[v]):
            if c in seen: continue
            seen.add(c)
            if c not in match or try_assign(match[c], seen):
                match[c] = v; return True
        return False
    for v in range(n):
        for _ in range(b):
            if not try_assign(v, set()): return False
    return True

ok = True
# 1. witness
L = [set(range(1,7)), {1,2,3,4,7,8}, {5,6,9,10,11,12}, {5,6,13,14,15,16}]
a, b = 6, 4
assert all(len(x) == a for x in L)
phi = [{1,2,5,6}, {3,4,7,8}, {9,10,11,12}, {13,14,15,16}]
valid = all(phi[i] <= L[i] and len(phi[i]) == b for i in range(4)) and all(not (phi[i] & phi[j]) for i in range(4) for j in range(i+1,4))
inter = [len(L[i] & L[j]) for i in range(4) for j in range(i+1,4)]
results = set()
for pick in (min, max):
    for perm_seed in range(3):
        w, out = colorsym(L, b, order_key=lambda S, s=perm_seed: (hash((tuple(sorted(S)), s)) % 97), pick=pick)
        results.add(tuple(w))
        okout = all(out[i] <= L[i] for i in range(4)) and all(not (out[i] & out[j]) for i in range(4) for j in range(i+1,4))
        ok &= okout
w0, out0 = colorsym(L, b)
line = valid and colorable(L, b) and results == {(5,3,4,4)} and w0 == [5,3,4,4]
ok &= line
print(f"witness n=4 a=6 b=4: lists ok, pairwise intersections {inter}, explicit coloring valid={valid}, "
      f"colorable={colorable(L,b)}, ColorSym loads={sorted(results)}, outputs (min-color) {[sorted(x) for x in out0]} {'OK' if line else 'FAIL'}")
# 2. family
for m in range(2, 13):
    P = set(range(1, 2*m+1)); T = set(range(2*m+1, 3*m+1)); Q2 = set(range(3*m+1, 4*m+1))
    Q3 = set(range(4*m+1, 6*m+1)); Q4 = set(range(6*m+1, 8*m+1))
    L = [P|T, P|Q2, T|Q3, T|Q4]; a, b = 3*m, 2*m
    assert all(len(x) == a for x in L)
    w, _ = colorsym(L, b)
    col = colorable(L, b)
    line = col and w[1] == m + m//2 and w[1] < b
    ok &= line
    print(f"family m={m}: a={a} b={b} colorable={col} loads={w} vertex2={w[1]} (m+floor(m/2)={m+m//2}) {'OK' if line else 'FAIL'}")
print("ALL_OK" if ok else "SOME_FAIL")
sys.exit(0 if ok else 1)
