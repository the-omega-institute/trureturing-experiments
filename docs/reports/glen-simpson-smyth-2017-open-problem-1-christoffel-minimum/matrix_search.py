"""Lengths reachable by compositions of A, B (L += 1 each) and h (L += 4 each) applied to ab (L = 3), using only
incidence matrices; report lengths n <= NMAX where the composition's factor count is below C(n).
The exact factor-count identities need every intermediate word to contain both letters (true for these images)."""
import sys, itertools
from math import gcd
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 506
MAXOPS = int(sys.argv[2]) if len(sys.argv) > 2 else 16
def qsum(a, b):
    s = 0
    while b: s += a // b; a, b = b, a % b
    return s
C = {n: 1 + min(qsum(p, n) for p in range(1, n) if gcd(p, n) == 1) for n in range(2, NMAX + 1)}
MA = ((1, 1), (0, 1)); MB = ((1, 0), (1, 1)); MH = ((2, 2), (1, 2))
def mul(X, Y): return tuple(tuple(sum(X[i][k] * Y[k][j] for k in range(2)) for j in range(2)) for i in range(2))
best = {}
# BFS over (matrix, L, word) with pruning by length
start = (((1, 0), (0, 1)), 3, '')
layer = [start]; seenM = {}
for depth in range(MAXOPS):
    nxt = []
    for M, L, w in layer:
        for name, G, inc in (('A', MA, 1), ('B', MB, 1), ('h', MH, 4)):
            if name == 'h' and w.count('h') >= 2: continue
            N = mul(M, G)                     # composition: existing outer word, new map applied first (innermost)
            n = sum(N[i][j] for i in range(2) for j in range(2))
            if n > NMAX: continue
            L2 = L + inc
            key = (N, w.count('h') + (name == 'h'))
            if seenM.get(key, 99) <= L2: continue
            seenM[key] = L2
            nxt.append((N, L2, w + name))
            if n not in best or L2 < best[n][0]: best[n] = (L2, w + name)
    layer = nxt
hits = sorted((n, best[n][0], C[n], best[n][1]) for n in best if n >= 2 and best[n][0] < C[n])
print(f"compositions of A, B and at most two h, at most {MAXOPS} maps, n <= {NMAX}: {len(hits)} counterexample lengths")
for h in hits[:20]: print("  n=%d L=%d C(n)=%d word=%s (outer→inner)" % h)
