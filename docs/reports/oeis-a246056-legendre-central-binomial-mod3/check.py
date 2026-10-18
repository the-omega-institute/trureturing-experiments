"""OEIS A246056 (P. D. Hanna): G.f. Sum_{n>=0} x^n/(1-2x)^(2n+1) [Sum_k C(n,k)^2 2^k x^k] [Sum_k C(n,k)^2 3^k x^k].
Conjecture: a(n) = 1 (mod 3) when n = 2*A005836(k), and a(n) = 0 (mod 3) otherwise (A005836: no digit 2 in base 3).

The program
1. computes a(n) from the g.f. exactly and compares with the entry's b-file;
2. checks the identity Sum_n x^n (1+x)^(-2n-1) Sum_k C(n,k)^2 (-x)^k = 1/sqrt(1+4x^2) = Sum_m (-1)^m C(2m,m) x^(2m)
   exactly (integer coefficients) to order N;
3. checks that the g.f. is congruent mod 3 termwise to that series (2 = -1, 3 = 0 mod 3), and the residues of a(n):
   (-1)^m C(2m,m) = 1 (mod 3) iff m has no base-3 digit 2 (Lucas), else 0.
Prints RESULT PASS or RESULT FAIL.
"""
import sys
from math import comb

N = int(sys.argv[1]) if len(sys.argv) > 1 else 300
BFILE = sys.argv[2] if len(sys.argv) > 2 else "b246056.txt"
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


def neg_binom_series(c, e, n):
    """(1 - c x)^(-e) to order n."""
    return [comb(e + k - 1, k) * c ** k for k in range(n + 1)]


def gf(c1, c2, c3, n):
    """Sum_m x^m (1 - c1 x)^(-2m-1) [Sum_k C(m,k)^2 c2^k x^k] [Sum_k C(m,k)^2 c3^k x^k]."""
    out = [0] * (n + 1)
    for m in range(n + 1):
        p2 = [comb(m, k) ** 2 * c2 ** k for k in range(m + 1)]
        p3 = [comb(m, k) ** 2 * c3 ** k for k in range(m + 1)]
        t = mul(mul(neg_binom_series(c1, 2 * m + 1, n - m), p2, n - m), p3, n - m)
        for i in range(n - m + 1):
            out[i + m] += t[i]
    return out


a = gf(2, 2, 3, N)
bf = {}
for line in open(BFILE):
    s = line.split()
    if len(s) == 2 and not line.startswith("#"):
        bf[int(s[0])] = int(s[1])
cmp = [k for k in range(0, N + 1) if k in bf]
check(all(a[k] == bf[k] for k in cmp), "a(0..%d) equal the b-file" % max(cmp))
red = gf(-1, -1, 0, N)                       # 2 -> -1, 3 -> 0, (1-2x) -> (1+x)
target = [0] * (N + 1)
for m in range(N // 2 + 1):
    target[2 * m] = (-1) ** m * comb(2 * m, m)
check(red == target, "Sum_n x^n (1+x)^(-2n-1) Sum_k C(n,k)^2 (-x)^k = Sum_m (-1)^m C(2m,m) x^(2m) to order %d" % N)
check(all((a[k] - red[k]) % 3 == 0 for k in range(N + 1)), "the g.f. is congruent mod 3 to the reduced series")


def no_two(m):
    while m:
        if m % 3 == 2:
            return False
        m //= 3
    return True


check(all(((-1) ** m * comb(2 * m, m)) % 3 == (1 if no_two(m) else 0) for m in range(N // 2 + 1)),
      "(-1)^m C(2m,m) = 1 (mod 3) iff m has no base-3 digit 2, else 0")
check(all((a[k] % 3) == (1 if (k % 2 == 0 and no_two(k // 2)) else 0) for k in range(N + 1)),
      "a(n) = 1 (mod 3) iff n = 2m with m in A005836, else 0, for n <= %d" % N)
print("RESULT " + ("PASS" if ok else "FAIL"))
sys.exit(0 if ok else 1)
