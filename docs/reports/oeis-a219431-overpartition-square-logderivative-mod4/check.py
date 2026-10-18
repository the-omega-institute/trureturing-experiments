"""OEIS A219431: logarithmic derivative of G(x) = sum_n pbar(n^2) x^n (A219430), pbar = overpartition counts (A015128).
Conjecture: a(2n+1) = 2 (mod 4), a(2n) = 0 (mod 4).

The program
1. computes pbar(m) for m <= M0 from the product prod_k (1 + q^k)/(1 - q^k) and checks pbar(m) = 2 d(m) (mod 4), m >= 1
   (the proof's step: prod_k (1 + 2 sum_j q^(kj)) = 1 + 2 sum_m d(m) q^m (mod 4));
2. computes pbar(m) for m <= N^2 from 1/phi(-q) (pbar(m) = -2 sum_{k>=1} (-1)^k pbar(m - k^2)), cross-checked
   against the product values for m <= M0, and G(x) = sum pbar(n^2) x^n;
3. computes a(n) = [x^n] x G'(x)/G(x) for n <= N, compares with the entry's b-file and checks a(n) = 2n (mod 4).
Prints RESULT PASS or RESULT FAIL.
"""
import sys

N = int(sys.argv[1]) if len(sys.argv) > 1 else 200
BFILE = sys.argv[2] if len(sys.argv) > 2 else "b219431.txt"
M0 = 3000
ok = True


def check(c, msg):
    global ok
    print(("ok   " if c else "FAIL ") + msg)
    ok = ok and c


# 1. product definition
P = [0] * (M0 + 1)
P[0] = 1
for k in range(1, M0 + 1):
    # multiply by (1 + q^k)
    for m in range(M0, k - 1, -1):
        P[m] += P[m - k]
    # divide by (1 - q^k)
    for m in range(k, M0 + 1):
        P[m] += P[m - k]
d = [0] * (M0 + 1)
for i in range(1, M0 + 1):
    for j in range(i, M0 + 1, i):
        d[j] += 1
check(all((P[m] - 2 * d[m]) % 4 == 0 for m in range(1, M0 + 1)), "pbar(m) = 2 d(m) (mod 4) for 1 <= m <= %d (product definition)" % M0)
# 2. theta recursion up to N^2
M = N * N
Q = [0] * (M + 1)
Q[0] = 1
for m in range(1, M + 1):
    s = 0
    k = 1
    while k * k <= m:
        s += (-1) ** k * Q[m - k * k]
        k += 1
    Q[m] = -2 * s
check(Q[: M0 + 1] == P, "the theta recursion equals the product for m <= %d" % M0)
G = [Q[n * n] for n in range(N + 1)]
check(all(G[n] % 4 == 2 for n in range(1, N + 1)), "pbar(n^2) = 2 (mod 4) for 1 <= n <= %d" % N)
# 3. logarithmic derivative: n G_n = sum_{k=1..n} a_k G_{n-k}
a = [0] * (N + 1)
for n in range(1, N + 1):
    a[n] = n * G[n] - sum(a[k] * G[n - k] for k in range(1, n))
bf = {}
for line in open(BFILE):
    s = line.split()
    if len(s) == 2 and not line.startswith("#"):
        bf[int(s[0])] = int(s[1])
cmp = [n for n in range(1, N + 1) if n in bf]
check(all(a[n] == bf[n] for n in cmp), "a(1..%d) equal the b-file" % max(cmp))
check(all((a[n] - 2 * n) % 4 == 0 for n in range(1, N + 1)), "a(n) = 2n (mod 4) for 1 <= n <= %d" % N)
print("RESULT " + ("PASS" if ok else "FAIL"))
sys.exit(0 if ok else 1)
