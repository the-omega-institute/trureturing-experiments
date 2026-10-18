"""OEIS A356508 (P. D. Hanna): 2 = Product_{n>=1} (1 + x^n A(x)) (1 + x^(n-1)/A(x)).
Conjecture: a(n) = 2 (mod 4) at n = 2k for every k > 0 with p(k) odd (A052002), otherwise a(n) = 0 (mod 4), n > 0.

The program
1. solves 1/A = 2/P - 1, P = Product_{n>=1} (1 + x^n A)(1 + x^n/A) (the factor n = 1 of the second product is
   1 + 1/A), by its division-free recursion, compares with the entry's b-file, and checks the defining product;
2. checks the proof steps: A = 1 (mod 2); each factor (1 + x^n A)(1 + x^n/A) = (1 + x^n)^2 (mod 4);
   1/A = 1 + 2 sum_{n>=1} p(n) x^(2n) (mod 4);
3. checks the residues of a(n) for 1 <= n <= N.
Prints RESULT PASS or RESULT FAIL.
"""
import sys

N = int(sys.argv[1]) if len(sys.argv) > 1 else 400
BFILE = sys.argv[2] if len(sys.argv) > 2 else "b356508.txt"
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


def inv_unit(p, n):
    r = [0] * (n + 1)
    r[0] = 1 if p[0] == 1 else None
    assert p[0] == 1
    for k in range(1, n + 1):
        r[k] = -sum(p[j] * r[k - j] for j in range(1, min(k, len(p) - 1) + 1))
    return r


def product_P(A, Ai, n):
    P = [1] + [0] * n
    for k in range(1, n + 1):
        f1 = [0] * (n + 1); f1[0] = 1
        f2 = [0] * (n + 1); f2[0] = 1
        for i in range(n + 1 - k):
            f1[i + k] += A[i]
            f2[i + k] += Ai[i]
        P = mul(P, mul(f1, f2, n), n)
    return P


A = [1] + [0] * N
for m in range(1, N + 1):
    # coefficient m of 1/A = 2/P - 1 uses A[0..m-1] only (each factor of P carries x^k, k >= 1)
    Ai = inv_unit(A[:m] + [0] * (N + 1 - m), m)
    P = product_P(A[:m] + [0] * (N + 1 - m), Ai, m)
    rhs = [2 * c for c in inv_unit(P, m)]
    rhs[0] -= 1
    target_inv = rhs                                   # 1/A to order m
    A[m] = inv_unit(target_inv, m)[m]
Ai = inv_unit(A, N)
P = product_P(A, Ai, N)
full = mul([1 + Ai[0]] + Ai[1:], P, N)                 # (1 + 1/A) * P
check(full == [2] + [0] * N, "2 = (1 + 1/A) prod_{n>=1} (1 + x^n A)(1 + x^n/A) to order %d" % N)
bf = {}
for line in open(BFILE):
    s = line.split()
    if len(s) == 2 and not line.startswith("#"):
        bf[int(s[0])] = int(s[1])
cmp = [k for k in range(0, N + 1) if k in bf]
check(all(A[k] == bf[k] for k in cmp), "a(0..%d) equal the b-file" % max(cmp))
check(A[0] == 1 and all(A[k] % 2 == 0 for k in range(1, N + 1)), "A = 1 (mod 2)")
ok_f = True
for k in range(1, 40):
    f1 = [0] * (N + 1); f1[0] = 1
    f2 = [0] * (N + 1); f2[0] = 1
    for i in range(N + 1 - k):
        f1[i + k] += A[i]
        f2[i + k] += Ai[i]
    g = mul(f1, f2, N)
    sq = [0] * (N + 1); sq[0] = 1
    if k <= N: sq[k] += 2
    if 2 * k <= N: sq[2 * k] += 1
    ok_f = ok_f and all((g[i] - sq[i]) % 4 == 0 for i in range(N + 1))
check(ok_f, "(1 + x^k A)(1 + x^k/A) = (1 + x^k)^2 (mod 4) for k < 40")
p = [1] + [0] * N
for k in range(1, N + 1):
    for m in range(k, N + 1):
        p[m] += p[m - k]
pred = [0] * (N + 1)
pred[0] = 1
for n in range(1, N // 2 + 1):
    pred[2 * n] = 2 * p[n]
check(all((Ai[i] - pred[i]) % 4 == 0 for i in range(N + 1)), "1/A = 1 + 2 sum p(n) x^(2n) (mod 4)")
check(all((A[m] % 4 == 2) if (m % 2 == 0 and p[m // 2] % 2 == 1) else (A[m] % 4 == 0) for m in range(1, N + 1)),
      "a(n) = 2 (mod 4) iff n = 2k with p(k) odd, else 0 (mod 4), 1 <= n <= %d" % N)
print("RESULT " + ("PASS" if ok else "FAIL"))
sys.exit(0 if ok else 1)
