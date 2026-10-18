"""OEIS A372527 (P. D. Hanna): A(A(x)^2) = x*A(x)/(1 - 2*A(x)), A = x + 2x^2 + ...
Conjecture: a(2^k + 1) = 2 (mod 4) for k >= 0, and a(n) = 0 (mod 4) for n > 0 not of the form 2^k + 1.

The program
1. solves the equation by its division-free integer recursion (a(n) enters the coefficient of x^(n+1) with
   multiplier 2 - 1 = 1) and compares with the entry's b-file;
2. checks the residues mod 4 for 2 <= n <= N, and records a(1) = 1;
3. checks the proof steps over F_2: x solves the equation mod 2, and C = (A - x)/2 satisfies
   C(x^2) = x*C + x^3 (mod 2), whose solution with c(1) = 0 is sum_{k>=0} x^(2^k + 1).
Prints RESULT PASS or RESULT FAIL.
"""
import sys

N = int(sys.argv[1]) if len(sys.argv) > 1 else 160
BFILE = sys.argv[2] if len(sys.argv) > 2 else "b372527.txt"
ok = True


def check(c, msg):
    global ok
    print(("ok   " if c else "FAIL ") + msg)
    ok = ok and c


def mul(p, q, n):
    r = [0] * (n + 1)
    for i, a in enumerate(p[: n + 1]):
        if a:
            for j in range(0, min(n - i, len(q) - 1) + 1):
                if q[j]:
                    r[i + j] += a * q[j]
    return r


def compose(a, b, n):
    r = [0] * (n + 1)
    for k in range(min(n, len(a) - 1), -1, -1):
        r = mul(r, b, n)
        r[0] += a[k]
    return r


def geom2(a, n):
    """1/(1 - 2a) for a with zero constant term (integer coefficients)."""
    r = [0] * (n + 1); r[0] = 1
    t = [2 * c for c in a[: n + 1]] + [0] * max(0, n + 1 - len(a))
    for k in range(1, n + 1):
        r[k] = sum(t[j] * r[k - j] for j in range(1, k + 1))
    return r


def residual(a, n):
    L = compose(a, mul(a, a, n), n)
    R = mul([0, 1], mul(a, geom2(a, n), n), n)
    return [L[i] - R[i] for i in range(n + 1)]


a = [0, 1] + [0] * N
for m in range(2, N + 1):
    a[m] = 0
    d = residual(a[: m + 1], m + 1)[m + 1]
    a[m] = -d
check(residual(a, N) == [0] * (N + 1), "A(A^2) = x A/(1 - 2A) to order %d" % N)
bf = {}
for line in open(BFILE):
    s = line.split()
    if len(s) == 2 and not line.startswith("#"):
        bf[int(s[0])] = int(s[1])
cmp = [n for n in range(1, N + 1) if n in bf]
check(all(a[n] == bf[n] for n in cmp), "a(1..%d) equal the b-file" % max(cmp))
pk = {2 ** k + 1 for k in range(0, 20)}
check(all((a[n] % 4 == 2) if n in pk else (a[n] % 4 == 0) for n in range(2, N + 1)),
      "a(n) = 2 (mod 4) iff n = 2^k + 1, else 0 (mod 4), for 2 <= n <= %d" % N)
check(a[1] == 1, "a(1) = 1 (the literal 'n > 0' clause fails only at n = 1)")
x = [0, 1] + [0] * N
check(all(c % 2 == 0 for c in residual(x, N)), "x solves the equation mod 2")
C = [(a[i] - x[i]) // 2 for i in range(N + 1)]
C2 = [0] * (N + 1)
for i, c in enumerate(C):
    if 2 * i <= N:
        C2[2 * i] = c
rhs = [0] * (N + 1)
for i in range(N):
    rhs[i + 1] += C[i]
rhs[3] += 1
check(all((C2[i] - rhs[i]) % 2 == 0 for i in range(N + 1)), "C(x^2) = x C + x^3 (mod 2)")
check(all((C[i] % 2) == (1 if i in pk and i >= 2 else 0) for i in range(N + 1)), "C = sum x^(2^k+1) (mod 2)")
print("RESULT " + ("PASS" if ok else "FAIL"))
sys.exit(0 if ok else 1)
