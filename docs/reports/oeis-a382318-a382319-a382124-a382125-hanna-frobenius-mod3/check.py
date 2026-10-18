"""OEIS A382318/A382319 and A382124/A382125 (P. D. Hanna): residues mod 3 through the Frobenius congruence.

A382318: A = x + (sum_{n>=1} A(x^n))^3.      A382319: B = x/(1-x) + sum_{n>=1} B(x^n)^3.
  Conjectures: a(3n) = b(n) (mod 3) for n >= 0; a(3n+1) = a(3n+2) = 0 (mod 3) for n >= 1.
  Proof steps checked: sum_{n>=1} A(x^n) = B, A = x + B^3, and B^3 = B(x^3) (mod 3).
A382124: C = exp(sum sigma(n)sigma(2n)/3 x^n/n).   A382125: D = exp(sum sigma(n)sigma(2n) x^n/n).
  Conjectures: c(n) = d(3n) (mod 3); d(3n+1) = d(3n+2) = 0 (mod 3).
  Proof steps checked: c_n = sigma(n)sigma(2n)/3 is an integer; the Dold congruences c_{pm} = c_m (mod p^(v_p(m)+1))
  for p <= P, m <= M; C has integer coefficients; D = C^3; C^3 = C(x^3) (mod 3).
Prints RESULT PASS or RESULT FAIL.
"""
import sys
from fractions import Fraction as Fr

N = int(sys.argv[1]) if len(sys.argv) > 1 else 300
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


def bfile(name):
    d = {}
    for line in open(name):
        s = line.split()
        if len(s) == 2 and not line.startswith("#"):
            d[int(s[0])] = int(s[1])
    return d


def lambert(a, n):
    """sum_{k>=1} a(x^k)."""
    r = [0] * (n + 1)
    for k in range(1, n + 1):
        for i in range(1, n // k + 1):
            r[i * k] += a[i]
    return r


# A382318 and A382319 by their own recursions (division-free)
a = [0] * (N + 1)
a[1] = 1
for m in range(2, N + 1):
    s = lambert(a, m)
    a[m] = mul(mul(s, s, m), s, m)[m] + (1 if m == 1 else 0)
b = [0] * (N + 1)
for m in range(1, N + 1):
    cube = 0
    for k in range(1, m + 1):
        if m % k == 0:
            bb = b[: m // k + 1]
            cube += mul(mul(bb, bb, m // k), bb, m // k)[m // k]
    b[m] = 1 + cube
b318, b319 = bfile("b382318.txt"), bfile("b382319.txt")
check(all(a[n] == b318[n] for n in range(1, N + 1) if n in b318), "A382318 a(1..%d) equal the b-file" % N)
check(all(b[n] == b319[n] for n in range(1, N + 1) if n in b319), "A382319 b(1..%d) equal the b-file" % N)
check(lambert(a, N) == b, "sum_{n>=1} A(x^n) = B to order %d" % N)
B3 = mul(mul(b, b, N), b, N)
check([a[i] - (1 if i == 1 else 0) for i in range(N + 1)] == B3, "A = x + B^3 to order %d" % N)
Bx3 = [0] * (N + 1)
for i in range(N // 3 + 1):
    Bx3[3 * i] = b[i]
check(all((B3[i] - Bx3[i]) % 3 == 0 for i in range(N + 1)), "B^3 = B(x^3) (mod 3) to order %d" % N)
check(all((a[3 * n] - b[n]) % 3 == 0 for n in range(0, N // 3 + 1)), "a(3n) = b(n) (mod 3), 3n <= %d" % N)
check(all(a[3 * n + 1] % 3 == 0 and a[3 * n + 2] % 3 == 0 for n in range(1, (N - 2) // 3 + 1)),
      "a(3n+1) = a(3n+2) = 0 (mod 3), 1 <= n, 3n+2 <= %d" % N)


# A382124 and A382125
def sigma(n):
    return sum(d for d in range(1, n + 1) if n % d == 0)


M = N
cs = [0] + [sigma(n) * sigma(2 * n) for n in range(1, 3 * M + 1)]
check(all(cs[n] % 3 == 0 for n in range(1, 3 * M + 1)), "sigma(n)sigma(2n) = 0 (mod 3) for n <= %d" % (3 * M))
c = [0] + [cs[n] // 3 for n in range(1, 3 * M + 1)]


def vp(m, p):
    v = 0
    while m % p == 0:
        m //= p
        v += 1
    return v


primes = [p for p in range(2, 40) if all(p % q for q in range(2, p))]
check(all((c[p * m] - c[m]) % p ** (vp(m, p) + 1) == 0 for p in primes for m in range(1, 3 * M // p + 1)),
      "Dold congruences c_{pm} = c_m (mod p^(v_p(m)+1)) for p < 40, pm <= %d" % (3 * M))


def expo(cc, n):
    """exp(sum cc_k x^k/k) via n e_n = sum_{k=1..n} cc_k e_{n-k}."""
    e = [Fr(1)] + [Fr(0)] * n
    for m in range(1, n + 1):
        e[m] = sum(cc[k] * e[m - k] for k in range(1, m + 1)) / m
    return e


C = expo(c, M)
D = expo(cs, 3 * M)
check(all(x.denominator == 1 for x in C), "C = A382124 has integer coefficients to order %d" % M)
check(all(x.denominator == 1 for x in D), "D = A382125 has integer coefficients to order %d" % (3 * M))
C = [int(x) for x in C]
D = [int(x) for x in D]
b124, b125 = bfile("b382124.txt"), bfile("b382125.txt")
check(all(C[n] == b124[n] for n in range(0, M + 1) if n in b124), "A382124 c(0..%d) equal the b-file" % M)
check(all(D[n] == b125[n] for n in range(0, 3 * M + 1) if n in b125), "A382125 d(0..) equal the b-file")
C3 = mul(mul(C, C, M), C, M)
check(C3 == D[: M + 1], "D = C^3 to order %d" % M)
check(all((D[3 * n] - C[n]) % 3 == 0 for n in range(0, M + 1)), "d(3n) = c(n) (mod 3), n <= %d" % M)
check(all(D[3 * n + 1] % 3 == 0 and D[3 * n + 2] % 3 == 0 for n in range(0, M)), "d(3n+1) = d(3n+2) = 0 (mod 3), n < %d" % M)
print("RESULT " + ("PASS" if ok else "FAIL"))
sys.exit(0 if ok else 1)
