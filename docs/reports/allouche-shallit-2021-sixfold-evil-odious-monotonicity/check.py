# Independent check of the sharp sixfold monotonicity of evil/odious representation counts
# (Allouche–Shallit arXiv:2112.13627v3, Conjecture 12) and of the prefix certificate route.
import sys
from math import comb
N = 1 << 13   # direct range
def t(n): return bin(n).count('1') & 1
E = [1 if t(n) == 0 else 0 for n in range(N)]
O = [1 - e for e in E]
def conv(a, b):
    r = [0]*N
    for i, x in enumerate(a):
        if x:
            for j in range(N - i):
                if b[j]: r[i+j] += x*b[j]
    return r
def power6(a):
    a2 = conv(a, a); a3 = conv(a2, a); return conv(a3, a3)
r6 = power6(E); s6 = power6(O)
# sharp thresholds
assert r6[36] > r6[37], (r6[36], r6[37])
assert all(r6[n] < r6[n+1] for n in range(37, N-1))
assert s6[4] == s6[5] == 0 and s6[6] == 1
assert all(s6[n] < s6[n+1] for n in range(5, N-1))
print('direct: r6 increasing on [37,%d), s6 on [5,%d); r6(36)=%d r6(37)=%d' % (N-1, N-1, r6[36], r6[37]))
# recurrences f(n) = sum_{k = n mod 2} p_k f((n-k)/2), f(0)=1
def poly_mul(a, b):
    r = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): r[i+j] += x*y
    return r
def ppow(p, e):
    r = [1]
    for _ in range(e): r = poly_mul(r, p)
    return r
P = {j: poly_mul(ppow([1, -1], j), ppow([1, 1], 5-j)) for j in range(0, 6)}
P[6] = ppow([1, -1], 6)
def table(p, M):
    f = [0]*M; f[0] = 1
    for n in range(1, M):
        f[n] = sum(pk * f[(n-k)//2] for k, pk in enumerate(p) if n-k >= 0 and (n-k) % 2 == 0)
    return f
M = 1 << 13
h = {j: table(P[j], M) for j in range(1, 6)}
c = table(P[6], M)
# h_0 = binom(n+4,4) coefficient of (1-z)^-5; check H_j generating identity by direct series for small range
T = [(-1)**t(n) for n in range(M)]
def series_conv(a, b, L):
    r = [0]*L
    for i in range(L):
        if a[i]:
            for j in range(L-i): r[i+j] += a[i]*b[j]
    return r
L = 600
ones = [1]*L
Tl = T[:L]
for j in range(1, 6):
    s = [1] + [0]*(L-1)
    for _ in range(j): s = series_conv(s, Tl, L)
    for _ in range(5-j): s = series_conv(s, ones, L)
    assert s == h[j][:L], j
s = [1] + [0]*(L-1)
for _ in range(6): s = series_conv(s, Tl, L)
assert s == c[:L]
print('recurrence tables match T^j/(1-z)^(5-j) and T^6 for n <', L)
def cget(a, n): return a[n] if n >= 0 else 0
for n in range(1, M):
    base = sum(comb(6, j) * cget(h[j], n) for j in range(1, 6)) + comb(n+4, 4) + cget(c, n) - cget(c, n-1)
    alt = sum((-1)**j * comb(6, j) * cget(h[j], n) for j in range(1, 6)) + comb(n+4, 4) + cget(c, n) - cget(c, n-1)
    if n < N:
        assert base == 64*(r6[n] - r6[n-1]), n
        assert alt == 64*(s6[n] - s6[n-1]), n
print('first-difference identities verified for 1 <= n <', N)
# state matrices and norms
def mats(p):
    d = len(p) - 1
    def coef(k): return p[k] if 0 <= k <= d else 0
    return [[[coef(2*a + b - i) for a in range(d)] for i in range(d)] for b in (0, 1)]
def mm(A, B): return [[sum(A[i][k]*B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]
def norm(A): return max(sum(abs(x) for x in row) for row in A)
for j in range(1, 6):
    print('j=%d max row-sum norm %d' % (j, max(norm(X) for X in mats(P[j]))))
C0, C1 = mats(P[6])
import itertools
ident = [[int(i == k) for k in range(6)] for i in range(6)]
norms = []
for L5 in range(6):
    best = 0
    for w in itertools.product((C0, C1), repeat=L5):
        A = ident
        for X in w: A = mm(X, A)
        best = max(best, norm(A))
    norms.append(best)
print('C word norms by length 0..5:', norms)
assert norms[5] < 16**5 and all(norms[l] <= 2 * 16**l for l in range(5))
# prefix certificate for several prefix lengths
def cert(bits):
    lo, hi = 1 << (bits-1), 1 << bits
    worst = None
    for p in range(lo, hi):
        Bp = sum(comb(6, j) * max(abs(cget(h[j], p-i)) for i in range(5)) for j in range(1, 6)) \
             + 4 * max(abs(cget(c, p-i)) for i in range(6))
        v = p**4 - 24*Bp
        if worst is None or v < worst[0]: worst = (v, p)
    return worst
for bits in range(8, 14):
    w = cert(bits)
    print('prefix bits %d: min p^4-24B_p = %d at p=%d' % (bits, w[0], w[1]))
print('ALL_OK')
