"""Evans-Francis-Jones, arXiv:2509.02874v1, Conjecture 4: the cut-off values alpha_n (roots of nP1) converge to 1/sqrt(5)
from above. nP1(alpha) = min_i K_{i,i+1} - max_i K_{i,i+2}, K = (I - alpha A_n)^{-1} - I for the path P_n.
Checks, from the definition only (no closed forms): every sign change of nP1 on a 60-point grid of (0, 1/2) for 6 <= n <= 30, 36, 40, the root
values, monotonicity in n, and the gap to 1/sqrt(5); and the source's Lemma 18 value nP~1(1/2)."""
from mpmath import mp, mpf, matrix, sqrt, inverse
mp.dps = 30
c = 1 / sqrt(5)

def p1(n, a):
    A = matrix(n, n)
    for i in range(n - 1): A[i, i + 1] = 1; A[i + 1, i] = 1
    M = matrix(n, n)
    for i in range(n):
        for j in range(n): M[i, j] = (1 if i == j else 0) - a * A[i, j]
    K = inverse(M)
    for i in range(n): K[i, i] -= 1
    return min(K[i, i + 1] for i in range(n - 1)) - max(K[i, i + 2] for i in range(n - 2))

def roots(n, grid=400):
    xs = [mpf(k) / (2 * grid) for k in range(1, grid)]   # (0, 1/2)
    vals = [p1(n, x) for x in xs]
    out = []
    for k in range(len(xs) - 1):
        if vals[k] == 0: out.append(xs[k])
        elif vals[k] * vals[k + 1] < 0:
            lo, hi, flo = xs[k], xs[k + 1], vals[k]
            for _ in range(90):
                mid = (lo + hi) / 2; fm = p1(n, mid)
                if fm * flo > 0: lo, flo = mid, fm
                else: hi = mid
            out.append((lo + hi) / 2)
    return out, vals[0], vals[-1]

prev = None; ok = True
for n in list(range(6, 31)) + [36, 40]:
    rs, v0, v1 = roots(n, grid=60)
    a = rs[0] if rs else None
    mono = (prev is None or (a is not None and a < prev))
    ok = ok and len(rs) == 1 and a > c and mono
    print(f"n={n:3d} roots_in_(0,1/2)={len(rs)} alpha_n={mp.nstr(a, 20)} alpha_n - 1/sqrt5={mp.nstr(a - c, 6)} decreasing={mono}", flush=True)
    prev = a
print("ALL: unique root, above 1/sqrt5, strictly decreasing:", ok)
