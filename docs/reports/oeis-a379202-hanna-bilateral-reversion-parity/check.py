"""OEIS A379202 (P. D. Hanna): 1/x = Sum_{n in Z} A(x)^n (A(x)^n + 2)^(n+1). Conjecture: a(n) is even for n > 1.

With S(y) = Sum_{n in Z} y^n (y^n + 2)^(n+1) (formula (5): A is the compositional inverse of B = 1/S), the program
1. computes S exactly: the terms n >= 0 are y^n (y^n + 2)^(n+1); the term n = -m is y^(m^2-2m) (1 + 2 y^m)^(1-m),
   checked against the literal Laurent expression with sympy for m = 1..6;
2. checks S = y^(-1) T with T an integer series, T = 1 (mod 2); B = y/T is integral with B = y (mod 2);
3. computes A as the compositional inverse of B, compares with the entry's b-file, checks the equation (1) in the form
   x * S(A(x)) = 1 to order N, and checks a(n) even for 1 < n <= N.
Prints RESULT PASS or RESULT FAIL.
"""
import sys

import sympy

N = int(sys.argv[1]) if len(sys.argv) > 1 else 120
BFILE = sys.argv[2] if len(sys.argv) > 2 else "b379202.txt"
ok = True


def check(c, msg):
    global ok
    print(("ok   " if c else "FAIL ") + msg)
    ok = ok and c


y = sympy.symbols("y")
check(all(sympy.simplify(y ** (-m) * (y ** (-m) + 2) ** (1 - m) - y ** (m * m - 2 * m) * (1 + 2 * y ** m) ** (1 - m)) == 0
          for m in range(1, 7)), "term n = -m equals y^(m^2-2m) (1+2y^m)^(1-m), m = 1..6")


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
    r[0] = 1
    for k in range(1, n + 1):
        r[k] = -sum(p[j] * r[k - j] for j in range(1, min(k, len(p) - 1) + 1))
    return r


def powr(p, e, n):
    r = [1] + [0] * n
    for _ in range(e):
        r = mul(r, p, n)
    return r


def T_series(F, n):
    """T = y * S(F) as a power series in x, for F = x + ... (F = y gives T(y)); order n."""
    out = [0] * (n + 1)
    # n = -1 term: F^(-1) -> y * F^(-1) = x/F
    xF = inv_unit([F[i + 1] for i in range(n + 1)] + [0], n)    # x/F with F/x unit
    for i in range(n + 1):
        out[i] += xF[i]
    k = 0
    while True:                                              # n >= 0 terms: x * F^k (F^k + 2)^(k+1), order k + 1
        if k + 1 > n:
            break
        Fk = powr(F, k, n)
        base = Fk[:]
        base[0] += 2
        t = mul(Fk, powr(base, k + 1, n), n)
        for i in range(n):
            out[i + 1] += t[i]
        k += 1
    m = 2
    while m * m - 2 * m + 1 <= n:                            # n = -m terms, m >= 2: x * F^(m^2-2m) (1 + 2F^m)^(1-m)
        Fm = powr(F, m, n)
        u = [2 * c for c in Fm]
        u[0] += 1
        t = mul(powr(F, m * m - 2 * m, n), powr(inv_unit(u, n), m - 1, n), n)
        for i in range(n - 0):
            if i + 1 <= n:
                out[i + 1] += t[i]
        m += 1
    return out


Y = [0, 1] + [0] * N
T = T_series(Y, N)
check(T[0] == 1 and all((T[i] % 2) == (1 if i == 0 else 0) for i in range(N + 1)), "T = y S(y) is an integer series = 1 (mod 2)")
B = [0] + inv_unit(T, N - 1)                                 # B = y / T
check(all(B[i] % 2 == (1 if i == 1 else 0) for i in range(N + 1)), "B = y/T = y (mod 2)")


def compose(a, b, n):
    r = [0] * (n + 1)
    for k in range(min(n, len(a) - 1), -1, -1):
        r = mul(r, b, n)
        r[0] += a[k]
    return r


# compositional inverse A of B
A = [0, 1] + [0] * (N - 1)
for k in range(2, N + 1):
    A[k] = 0
    c = compose(B, A, k)
    A[k] = -c[k]
check(compose(B, A, N) == [0, 1] + [0] * (N - 1), "B(A(x)) = x to order %d" % N)
bf = {}
for line in open(BFILE):
    s = line.split()
    if len(s) == 2 and not line.startswith("#"):
        bf[int(s[0])] = int(s[1])
cmp = [k for k in range(1, N + 1) if k in bf]
check(all(A[k] == bf[k] for k in cmp), "a(1..%d) equal the b-file" % max(cmp))
TA = T_series(A, N - 1)                                      # x S(A) as a series in x
check(TA == [1] + [0] * (N - 1), "x * S(A(x)) = 1 to order %d, i.e. equation (1)" % (N - 1))
check(all(A[k] % 2 == 0 for k in range(2, N + 1)), "a(n) even for 1 < n <= %d" % N)
print("RESULT " + ("PASS" if ok else "FAIL"))
sys.exit(0 if ok else 1)
