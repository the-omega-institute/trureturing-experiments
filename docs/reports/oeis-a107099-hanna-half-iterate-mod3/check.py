"""OEIS A107099 (P. D. Hanna): A(A(x)) = x + 4x^3, A(x) = sum_{n>=0} a(n) x^(2n+1).
Conjecture: [x^n] A(x) = 0 (mod 3) except at n = 3^k.

The program
1. writes A = x + 2B and solves B + B(x + 2B) ... i.e. B(x) + B(x + 2B(x)) = 2x^3 by an integer, division-free
   recursion (B(x + 2B) - B(x) = sum_{j>=1} 2^j B^j * B^{(j)}/j!, and b_m enters only the linear term);
2. checks A(A(x)) = x + 4x^3 to order N and compares a(0..) with the b-file b107099.txt;
3. checks that [x^n]A = 0 (mod 3) unless n = 3^k, and that [x^(3^k)]A = binom(1/2, k) (mod 3), i.e. that A mod 3
   is the additive series whose Frobenius symbol is (1 + F)^(1/2).
Prints RESULT PASS or RESULT FAIL.
"""
import sys
from fractions import Fraction as Fr
from math import comb

N = int(sys.argv[1]) if len(sys.argv) > 1 else 243
BFILE = sys.argv[2] if len(sys.argv) > 2 else "b107099.txt"
ok = True


def check(c, msg):
    global ok
    print(("ok   " if c else "FAIL ") + msg)
    ok = ok and c


def mul(p, q, n):
    r = [0] * (n + 1)
    for i, a in enumerate(p[: n + 1]):
        if a:
            for j in range(0, n + 1 - i):
                if q[j]:
                    r[i + j] += a * q[j]
    return r


def compose(a, b, n):
    """a(b(x)) for b with zero constant term (Horner)."""
    r = [0] * (n + 1)
    for k in range(n, -1, -1):
        r = mul(r, b, n)
        r[0] += a[k]
    return r


# B: b_m for odd m >= 3. E(B) := B(x + 2B(x)) + B(x) - 2x^3 = 0.
B = [0] * (N + 1)
for m in range(3, N + 1, 2):
    B[m] = 0
    x2B = [0] * (N + 1); x2B[1] = 1
    for i in range(N + 1):
        x2B[i] += 2 * B[i]
    E = compose(B, x2B, m)
    E = [E[i] + B[i] for i in range(m + 1)]
    E[3] -= 2
    # b_m enters E at degree m with coefficient 2 (B(x) and the leading term of B(x + 2B)); E[m] = 2 b_m + rest
    rest = E[m]
    assert rest % 2 == 0, m
    B[m] = -rest // 2
A = [0] * (N + 1); A[1] = 1
for i in range(N + 1):
    A[i] += 2 * B[i]
AA = compose(A, A, N)
tgt = [0] * (N + 1); tgt[1] = 1; tgt[3] = 4
check(AA == tgt, "A(A(x)) = x + 4x^3 to order %d (integer coefficients)" % N)
bf = {}
for line in open(BFILE):
    s = line.split()
    if len(s) == 2 and not line.startswith("#"):
        bf[int(s[0])] = int(s[1])
cmp = [n for n in range(0, (N - 1) // 2 + 1) if n in bf]
check(all(A[2 * n + 1] == bf[n] for n in cmp), "a(0..%d) equal the b-file" % max(cmp))
check(all(A[i] == 0 for i in range(0, N + 1, 2)), "A is odd")
pow3 = {3 ** k for k in range(0, 12)}
check(all(A[i] % 3 == 0 for i in range(1, N + 1) if i not in pow3), "[x^n]A = 0 (mod 3) for n <= %d not a power of 3" % N)


def binhalf_mod3(k):
    v = Fr(1)
    for j in range(k):
        v *= Fr(1, 2) - j
    for j in range(1, k + 1):
        v /= j
    return v.numerator * pow(v.denominator, -1, 3) % 3


ks = [k for k in range(0, 12) if 3 ** k <= N]
check(all(A[3 ** k] % 3 == binhalf_mod3(k) for k in ks), "[x^(3^k)]A = binom(1/2,k) (mod 3) for 3^k <= %d" % N)
print("residues at 3^k:", [(3 ** k, A[3 ** k] % 3) for k in ks])
print("RESULT " + ("PASS" if ok else "FAIL"))
sys.exit(0 if ok else 1)
