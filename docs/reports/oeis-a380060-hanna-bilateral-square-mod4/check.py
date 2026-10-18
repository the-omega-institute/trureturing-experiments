"""OEIS A380060 (P. D. Hanna): A(x) = Sum_{n in Z} x^(n^2) (A(x) - x^n)^n / (1 - x^n A(x))^n.
Conjecture: for n > 0, a(n) = 2 (mod 4) iff n is a square, and 4 | a(n) otherwise.

The program
1. checks with exact rational functions that the n and -n summands coincide (Laurent identity, n = 1..6);
2. solves A = 1 + 2 sum_{m>=1} x^(m^2) ((A - x^m)/(1 - x^m A))^m by its division-free recursion and compares
   with the entry's b-file;
3. checks the residues mod 4 and the proof step: with A = 1 + 2C, ((A - x^m)/(1 - x^m A)) = 1 (mod 2), so
   C = sum_{m>=1} x^(m^2) (mod 2).
Prints RESULT PASS or RESULT FAIL.
"""
import sys

import sympy

N = int(sys.argv[1]) if len(sys.argv) > 1 else 400
BFILE = sys.argv[2] if len(sys.argv) > 2 else "b380060.txt"
ok = True


def check(c, msg):
    global ok
    print(("ok   " if c else "FAIL ") + msg)
    ok = ok and c


x, a = sympy.symbols("x a")
same = all(sympy.simplify(x ** (n * n) * (a - x ** n) ** n / (1 - x ** n * a) ** n
                          - x ** (n * n) * (a - x ** (-n)) ** (-n) / (1 - x ** (-n) * a) ** (-n)) == 0 for n in range(1, 7))
check(same, "the n and -n summands coincide as rational functions (n = 1..6)")


def mul(p, q, n):
    r = [0] * (n + 1)
    for i, u in enumerate(p[: n + 1]):
        if u:
            for j in range(0, min(n - i, len(q) - 1) + 1):
                if q[j]:
                    r[i + j] += u * q[j]
    return r


def inv_unit(p, n):
    """1/p for an integer series with p[0] = 1."""
    r = [0] * (n + 1)
    r[0] = 1
    for k in range(1, n + 1):
        r[k] = -sum(p[j] * r[k - j] for j in range(1, min(k, len(p) - 1) + 1))
    return r


def rhs(A, n):
    out = [0] * (n + 1)
    out[0] = 1
    m = 1
    while m * m <= n:
        num = A[:]; num[m] -= 1                          # A - x^m
        den = [0] * (n + 1); den[0] = 1
        for i in range(n + 1 - m):
            den[i + m] -= A[i]                           # 1 - x^m A
        ratio = mul(num, inv_unit(den, n), n)
        pw = [1] + [0] * n
        for _ in range(m):
            pw = mul(pw, ratio, n - m * m)
        for i in range(n - m * m + 1):
            out[i + m * m] += 2 * pw[i]
        m += 1
    return out


A = [0] * (N + 1)
A[0] = 1
for k in range(1, N + 1):
    A[k] = rhs(A[: k] + [0] * (N + 1 - k), k)[k]          # coefficient k uses A[0..k-1] only
check(rhs(A, N) == A, "A = 1 + 2 sum x^(m^2) ((A - x^m)/(1 - x^m A))^m to order %d" % N)
bf = {}
for line in open(BFILE):
    s = line.split()
    if len(s) == 2 and not line.startswith("#"):
        bf[int(s[0])] = int(s[1])
cmp = [k for k in range(0, N + 1) if k in bf]
check(all(A[k] == bf[k] for k in cmp), "a(0..%d) equal the b-file" % max(cmp))
sq = {k * k for k in range(1, 100)}
check(all((A[k] % 4 == 2) if k in sq else (A[k] % 4 == 0) for k in range(1, N + 1)),
      "for 1 <= n <= %d: a(n) = 2 (mod 4) iff n is a square, else 0 (mod 4)" % N)
C = [(A[k] - (1 if k == 0 else 0)) // 2 for k in range(N + 1)]
check(all(C[k] % 2 == (1 if k in sq else 0) for k in range(N + 1)), "C = sum x^(m^2) (mod 2)")
print("RESULT " + ("PASS" if ok else "FAIL"))
sys.exit(0 if ok else 1)
