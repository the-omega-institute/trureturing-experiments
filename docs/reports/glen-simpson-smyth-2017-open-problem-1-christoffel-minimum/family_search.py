"""Search the proved family x = P h Q(ab) (P, Q words over the elementary morphisms A, B; L(x) = 7 + |P| + |Q|)
for lengths n with 7 + |P| + |Q| < C(n) = 1 + min over p coprime to n of the partial-quotient sum of p/n.
Also directly verify the smallest hits (Lyndon, exact factor count, Christoffel minimum by direct counting)."""
import itertools, sys
from math import gcd
exec(open(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), 'verify506.py')).read().split("def main")[0])

MA = ((1, 1), (0, 1)); MB = ((1, 0), (1, 1)); MH = ((2, 2), (1, 2))
def mul(X, Y): return tuple(tuple(sum(X[i][k] * Y[k][j] for k in range(2)) for j in range(2)) for i in range(2))
def length(M, vec=(1, 1)):  # |f(ab)| where f has incidence M (rows: output letters, cols: source letters); source ab has counts (1,1)
    return sum(M[i][j] * vec[j] for i in range(2) for j in range(2))
Cmin = {}
def C(n):
    if n not in Cmin: Cmin[n] = 1 + min(qsum(p, n) for p in range(1, n) if gcd(p, n) == 1)
    return Cmin[n]
best = {}
MAXOPS = 13
mats = {0: [((), ((1, 0), (0, 1)))]}
for r in range(1, MAXOPS + 1):
    mats[r] = [(w + (c,), mul(M, MA if c == 'A' else MB)) for (w, M) in mats[r - 1] for c in 'AB']
hits = []
for r in range(0, MAXOPS + 1):
    for s in range(0, MAXOPS + 1 - r):
        for (P, MP) in mats[r]:
            for (Q, MQ) in mats[s]:
                M = mul(mul(MP, MH), MQ); n = length(M)
                if n > 4000: continue
                Lx = 7 + r + s
                if Lx < C(n): hits.append((n, Lx, C(n), ''.join(P), ''.join(Q)))
hits.sort()
print(f"family x = P h Q(ab), |P| + |Q| <= {MAXOPS}, n <= 4000: {len(hits)} hits with L(x) < C(n)")
seen = set()
for h in hits[:15]:
    if h[0] in seen: continue
    seen.add(h[0]); print("  n=%d L(x)=%d C(n)=%d P=%s Q=%s" % h)
# direct verification of the smallest hit
n, Lx, Cn, P, Q = hits[0]
app = lambda m, x: "".join(m[c] for c in x)
Amap = {"a": "a", "b": "ab"}; Bmap = {"a": "ab", "b": "b"}; Hmap = {"a": "aab", "b": "aabb"}
x = "ab"
for c in reversed(Q): x = app(Amap if c == 'A' else Bmap, x)
x = app(Hmap, x)
for c in reversed(P): x = app(Amap if c == 'A' else Bmap, x)
direct_C = min(L_fast(christoffel(p, n)) for p in range(1, n) if gcd(p, n) == 1)
print(f"smallest hit n={n}: |x|={len(x)} Lyndon={is_lyndon(x)} L(x) brute={L_brute(x)[0]} (claimed {Lx}); "
      f"Christoffel minimum by direct counting = {direct_C} (formula {Cn})")
ok = len(x) == n and is_lyndon(x) and L_brute(x)[0] == Lx and direct_C == Cn and Lx < Cn
print("SMALLEST", n, x if n < 400 else x[:60] + "...")
print("RESULT", "PASS" if ok else "FAIL"); sys.exit(0 if ok else 1)
