"""OEIS A246423 (P. D. Hanna): G.f. Sum_{n>=0} x^n/(1-3x)^(2n+1) [Sum_k C(n,k)^2 2^k x^k] [Sum_k C(n,k)^2 3^k x^k].
Conjecture: a(n) = 2 (mod 4) iff n = 2^k for k >= 2, and a(n) = 0 (mod 4) elsewhere except a(0) = 1.

The program
1. computes a(n) from the g.f. exactly and compares with the entry's b-file;
2. checks the exact identities R = Sum_n x^n (1+x)^(-2n-1) Q_n(-x) = 1/sqrt(1+4x^2) and
   T = Sum_n n x^n (1+x)^(-2n-1) Q_n(-x) = x (1-2x) (1+4x^2)^(-3/2), Q_n(t) = Sum_k C(n,k)^2 t^k;
3. checks the reduction a = R + 2x T' (mod 4) with T' = Sum_n n^2 x^n (1+x)^(-2n-1) Q_n(-x), T' = T (mod 2), and the
   residues of a(n).
Prints RESULT PASS or RESULT FAIL.
"""
import sys
from fractions import Fraction as Fr
from math import comb

N = int(sys.argv[1]) if len(sys.argv) > 1 else 300
BFILE = sys.argv[2] if len(sys.argv) > 2 else "b246423.txt"
ok = True


def check(c, msg):
    global ok
    print(("ok   " if c else "FAIL ") + msg)
    ok = ok and c


def mul(p, q, n):
    r = [0] * (n + 1)
    for i, u in enumerate(p[: n + 1]):
        if u:
            for j in range(0, min(n - i, len(q) - 1) + 1):
                if q[j]:
                    r[i + j] += u * q[j]
    return r


def nb(c, e, n):
    return [comb(e + k - 1, k) * c ** k for k in range(n + 1)]


def gf(c1, c2, c3, n, weight=lambda m: 1):
    out = [0] * (n + 1)
    for m in range(n + 1):
        p2 = [comb(m, k) ** 2 * c2 ** k for k in range(m + 1)]
        p3 = [comb(m, k) ** 2 * c3 ** k for k in range(m + 1)]
        t = mul(mul(nb(c1, 2 * m + 1, n - m), p2, n - m), p3, n - m)
        w = weight(m)
        for i in range(n - m + 1):
            out[i + m] += w * t[i]
    return out


a = gf(3, 2, 3, N)
bf = {}
for line in open(BFILE):
    s = line.split()
    if len(s) == 2 and not line.startswith("#"):
        bf[int(s[0])] = int(s[1])
cmp = [k for k in range(0, N + 1) if k in bf]
check(all(a[k] == bf[k] for k in cmp), "a(0..%d) equal the b-file" % max(cmp))
def series(weight, n):
    out = [0] * (n + 1)
    for m in range(n + 1):
        q = [comb(m, k) ** 2 * (-1) ** k for k in range(m + 1)]
        t = mul(nb(-1, 2 * m + 1, n - m), q, n - m)
        w = weight(m)
        for i in range(n - m + 1):
            out[i + m] += w * t[i]
    return out


R = series(lambda m: 1, N)
T = series(lambda m: m, N)
T2 = series(lambda m: m * m, N)
targetR = [0] * (N + 1)
for m in range(N // 2 + 1):
    targetR[2 * m] = (-1) ** m * comb(2 * m, m)
check(R == targetR, "R = Sum (-1)^m C(2m,m) x^(2m) = 1/sqrt(1+4x^2) to order %d" % N)
# (1+4x^2)^(-3/2) = sum binom(-3/2,k) 4^k x^(2k)
g = [Fr(0)] * (N + 1)
for k in range(N // 2 + 1):
    b = Fr(1)
    for j in range(k):
        b *= Fr(-3, 2) - j
    for j in range(1, k + 1):
        b /= j
    g[2 * k] = b * 4 ** k
tt = mul([0, 1, -2], g, N)
check([Fr(v) for v in T] == tt, "T = x(1-2x)(1+4x^2)^(-3/2) to order %d" % N)
check(all((T2[i] - T[i]) % 2 == 0 for i in range(N + 1)), "T' = T (mod 2)")
red = [R[i] + (2 * T2[i - 1] if i >= 1 else 0) for i in range(N + 1)]
check(all((a[i] - red[i]) % 4 == 0 for i in range(N + 1)), "a = R + 2x T' (mod 4)")
pw = {2 ** k for k in range(2, 20)}
check(a[0] == 1 and all((a[n] % 4 == 2) if n in pw else (a[n] % 4 == 0) for n in range(1, N + 1)),
      "a(n) = 2 (mod 4) iff n = 2^k, k >= 2, else 0 (mod 4), n >= 1; a(0) = 1")
print("RESULT " + ("PASS" if ok else "FAIL"))
sys.exit(0 if ok else 1)
