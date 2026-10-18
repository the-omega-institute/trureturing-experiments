"""OEIS A363558, A363559, A363561, A363569 (P. D. Hanna): coefficients a(m) of
A(x) = Sum_{n=-oo..+oo} x^n (c + x^n)^(2n) for c = 2, 3, sqrt(2), i. Conjectures, for each c:
  a(2^n) = c^(2^(n+1)) for n > 0, and a(p) = p c^(p-1) + c^(2p) for primes p > 3
(stated in the entries as 4^(2^n), 9^(2^n), 2^(2^n), 1 and p 2^(p-1) + 4^p, p 3^(p-1) + 9^p,
p 2^((p-1)/2) + 2^p, p (-1)^((p-1)/2) - 1).

Proof checked here. The term n >= 1 contributes C(2n, m/n - 1) c^(2n - m/n + 1) to x^m when n | m; the term
n = -k (k >= 1) equals x^(2k^2 - k) (1 + c x^k)^(-2k) and contributes C(-2k, m/k - 2k + 1) c^(m/k - 2k + 1) when k | m.
- m = p prime > 3: only n = p (c^(2p)) and k = 1 (C(-2, p-1) c^(p-1) = p c^(p-1)) remain.
- m = 2^N: the positive term at n = 2^i (i < N) and the negative term at k = 2^(N-1-i) have the same power of c and
  opposite coefficients (C(-a, b) = (-1)^b C(a+b-1, b) with b odd), so only n = 2^N (c^(2^(N+1))) remains.
The script compares the divisor formula with a direct series expansion (exact, sympy, symbolic c), with the published
data, and checks both statements for primes p <= P and 2^N <= 2^NMAX, as polynomials in c.
"""
import sys
from sympy import symbols, binomial, divisors, Rational, expand, sqrt, I, primerange, series, Poly
P = int(sys.argv[1]) if len(sys.argv) > 1 else 300
NMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 12
c, x = symbols('c x')
def gbinom(a, b):
    r = Rational(1)
    for t in range(b): r = r * (a - t) / (t + 1)
    return r
def coef(m, cc=c):
    s = 0
    for n in divisors(m):
        j = m // n - 1
        if 0 <= j <= 2 * n: s += binomial(2 * n, j) * cc ** (2 * n - j)
    for k in divisors(m):
        i = m // k - 2 * k + 1
        if i >= 0: s += gbinom(-2 * k, i) * cc ** i
    return expand(s)
# direct expansion of the bilateral sum up to x^M
M = 14
direct = 0
for n in range(0, M + 1):
    direct += x ** n * (c + x ** n) ** (2 * n)
for k in range(1, M + 1):
    if 2 * k * k - k <= M:
        direct += x ** (2 * k * k - k) * series((1 + c * x ** k) ** (-2 * k), x, 0, M + 1).removeO()
direct = expand(direct)
ok_series = all(expand(direct.coeff(x, m) - (coef(m) if m > 0 else 1)) == 0 for m in range(0, M + 1))
DATA = {2: [1, 5, 16, 77, 256, 1104, 4121, 16832, 65536, 264688, 1048617, 4205568, 16779008],
        I: [1, 0, 1, -3, 1, 4, -4, -8, 1, 23, -8, -12, -27]}
ok_data = all(expand(coef(m, cc)) == DATA[cc][m] for cc in DATA for m in range(1, len(DATA[cc])))
ok_pow = all(expand(coef(2 ** N) - c ** (2 ** (N + 1))) == 0 for N in range(1, NMAX + 1))
ok_prime = all(expand(coef(p) - (p * c ** (p - 1) + c ** (2 * p))) == 0 for p in primerange(5, P + 1))
inst = all(expand(coef(p, cc) - v) == 0 for p in primerange(5, 60) for cc, v in
           ((2, p * 2 ** (p - 1) + 4 ** p), (3, p * 3 ** (p - 1) + 9 ** p), (sqrt(2), p * 2 ** ((p - 1) // 2) + 2 ** p), (I, p * (-1) ** ((p - 1) // 2) - 1)))
print(f"divisor formula = direct expansion (symbolic c), m<={M}: {ok_series}; data (c = 2, i): {ok_data}")
print(f"a(2^N) = c^(2^(N+1)), 1<=N<={NMAX}: {ok_pow}; a(p) = p c^(p-1) + c^(2p), 5<=p<={P}: {ok_prime}; instances c=2,3,sqrt2,i: {inst}")
good = ok_series and ok_data and ok_pow and ok_prime and inst
print("RESULT", "PASS" if good else "FAIL")
sys.exit(0 if good else 1)
