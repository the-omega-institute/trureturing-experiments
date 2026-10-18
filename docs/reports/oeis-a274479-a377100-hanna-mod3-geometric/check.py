"""OEIS A274479 and A377100 (P. D. Hanna): a(n) = 1 (mod 3) for n >= 1.

A274479: A(x)^2 = A(x^2/(1-2x-4x^2)), A = x + ...
A377100: A(x) = A(x^3)/A(x^2) + A(x)^2, A = x + ...

For each entry the program
1. solves the functional equation coefficientwise over Q to order N and compares with the entry's example terms;
2. checks that every coefficient is an integer and = 1 (mod 3);
3. checks the proof steps: for A274479 the reversion R satisfies R(y^2)(1 - 2R - 4R^2) = R^2 and every division by 2
   in its recursion is exact; for both entries x/(1-x) satisfies the equation modulo 3.
Prints RESULT PASS or RESULT FAIL.
"""
from fractions import Fraction as Fr
import sys

N = int(sys.argv[1]) if len(sys.argv) > 1 else 200
ok = True


def check(cond, msg):
    global ok
    print(("ok   " if cond else "FAIL ") + msg)
    ok = ok and cond


def mul(p, q, n):
    r = [0] * (n + 1)
    for i, a in enumerate(p[: n + 1]):
        if a:
            for j, b in enumerate(q[: n + 1 - i]):
                if b:
                    r[i + j] += a * b
    return r


def inv(p, n):
    """Inverse of a series with constant term 1 (or a unit in the coefficient ring)."""
    r = [0] * (n + 1)
    r[0] = 1 / p[0] if not isinstance(p[0], int) or p[0] not in (1, -1) else p[0]
    for k in range(1, n + 1):
        r[k] = -sum(p[j] * r[k - j] for j in range(1, min(k, len(p) - 1) + 1)) * r[0]
    return r


def compose(a, b, n):
    """a(b(x)) for b with zero constant term."""
    r = [0] * (n + 1)
    pw = [1] + [0] * n
    for k in range(0, n + 1):
        if k > 0:
            pw = mul(pw, b, n)
        if k < len(a) and a[k]:
            for i in range(n + 1):
                r[i] += a[k] * pw[i]
    return r


def is_int(v):
    return isinstance(v, int) or v.denominator == 1


# ---------- A274479 ----------
phi = [0, 0] + inv([1, -2, -4], N)[: N - 1]           # x^2/(1-2x-4x^2)
phi = [Fr(c) for c in phi]
pw = {1: phi[:]}
for k in range(2, N // 2 + 2):
    pw[k] = mul(pw[k - 1], phi, N)
a = [Fr(0)] * (N + 1)
a[1] = Fr(1)
for m in range(2, N + 1):
    n = m + 1
    if n > N:
        break
    known = sum(a[i] * a[n - i] for i in range(1, n) if i != m and n - i != m)
    rhs = sum(a[k] * pw[k][n] for k in range(1, n // 2 + 1))
    a[m] = (rhs - known) / 2
M = N - 1                                              # a[1..M] determined
A1 = a[: M + 1]
ex1 = [1, 1, 4, 10, 34, 106, 361, 1219, 4252, 14932, 53263, 191533]
check(A1[1:13] == ex1, "A274479 first 12 terms equal the entry's example")
check(all(is_int(A1[n]) for n in range(1, M + 1)), "A274479 a(1..%d) are integers" % M)
check(all(int(A1[n]) % 3 == 1 for n in range(1, M + 1)), "A274479 a(1..%d) = 1 (mod 3)" % M)
lhs = mul(A1, A1, M)
rhs = compose(A1, phi, M)
check(lhs[: M + 1] == rhs[: M + 1], "A274479 A^2 = A(x^2/(1-2x-4x^2)) to order %d" % M)
# reversion recursion with integer arithmetic: R(y^2)(1-2R-4R^2) = R^2, solve [y^(n+1)] for r_n
R = [0, 1] + [0] * (M - 1)
exact = True
for nn in range(2, M):
    Rs = R[:]                                          # r_nn currently 0
    R2 = [0] * (M + 1)
    for i, c in enumerate(Rs):
        if 2 * i <= M:
            R2[2 * i] = c
    t = [1] + [0] * M
    t = [t[i] - 2 * Rs[i] - 4 * mul(Rs, Rs, M)[i] for i in range(M + 1)]
    e = [mul(R2, t, M)[i] - mul(Rs, Rs, M)[i] for i in range(M + 1)]
    # r_nn enters only R^2 at degree nn+1 with coefficient 2 (the other occurrences have degree >= nn+2)
    num = e[nn + 1]
    if num % 2:
        exact = False
    R[nn] = num // 2
check(exact, "A274479 every division by 2 in the reversion recursion is exact (n < %d)" % M)
comp = compose(A1, R, M - 1)
check(comp[: M] == [0, 1] + [0] * (M - 2), "A274479 R is the compositional inverse of A to order %d" % (M - 1))
# F = x/(1-x) solves the equation over F_3
F = [0] + [1] * M
lhsF = [c % 3 for c in mul(F, F, M)]
rhsF = [int(c) % 3 for c in compose(F, [int(c) for c in phi], M)]
check(lhsF == rhsF, "A274479 x/(1-x) satisfies the equation mod 3 to order %d" % M)

# ---------- A377100 ----------
b = [0, 1] + [0] * (N - 1)
for n in range(2, N + 1):
    # A(x^3)/A(x^2) = x * T(x^3)/T(x^2), T = A/x; coefficient of x^n uses b_k with k <= (n+1)/2
    T = [b[k + 1] for k in range(N)]
    T3 = [0] * (N + 1)
    T2 = [0] * (N + 1)
    for k, c in enumerate(T):
        if 3 * k <= N:
            T3[3 * k] = c
        if 2 * k <= N:
            T2[2 * k] = c
    q = mul(T3, inv(T2, N), N)
    sq = sum(b[i] * b[n - i] for i in range(1, n))
    b[n] = q[n - 1] + sq
ex2 = [1, 1, 1, 4, 10, 28, 82, 256, 808, 2617, 8566, 28447]
check(b[1:13] == ex2, "A377100 first 12 terms equal the entry's example")
check(all(isinstance(b[n], int) for n in range(1, N + 1)), "A377100 a(1..%d) are integers (division-free recursion)" % N)
check(all(b[n] % 3 == 1 for n in range(1, N + 1)), "A377100 a(1..%d) = 1 (mod 3)" % N)
lhs = [b[i] * 1 for i in range(N + 1)]
T = [b[k + 1] for k in range(N)]
T3 = [0] * (N + 1); T2 = [0] * (N + 1)
for k, c in enumerate(T):
    if 3 * k <= N: T3[3 * k] = c
    if 2 * k <= N: T2[2 * k] = c
quot = [0] + mul(T3, inv(T2, N), N)[:N]
rhs = [quot[i] + mul(b, b, N)[i] for i in range(N + 1)]
check(lhs == rhs, "A377100 A = A(x^3)/A(x^2) + A^2 to order %d" % N)
# x/(1-x) mod 3: (x + 2x^2)/(1-x)^2 = x/(1-x) mod 3, i.e. x(1-x^2)/(1-x^3) + x^2/(1-x)^2 = x/(1-x) mod 3
Fq = [0] + [((1 if (i - 1) % 3 == 0 else 0) - (1 if (i - 1) % 3 == 2 else 0)) for i in range(1, N + 1)]  # x(1-x^2)/(1-x^3)
F2 = [0, 0] + [k - 1 for k in range(2, N + 1)]          # x^2/(1-x)^2
check(all((Fq[i] + F2[i] - (1 if i >= 1 else 0)) % 3 == 0 for i in range(N + 1)),
      "A377100 x/(1-x) satisfies the equation mod 3 to order %d" % N)

print("RESULT " + ("PASS" if ok else "FAIL"))
sys.exit(0 if ok else 1)
