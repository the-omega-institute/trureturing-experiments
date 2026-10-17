"""OEIS A386655, A386656, A386657, A386658 (P. D. Hanna): e.g.f. Sum_{n>=0} (q^n x + LambertW(x))^n / n! for
q = 2, 3, 4, 5. Conjecture: for n >= 6, a(n) mod 6 repeats with period 6 (patterns [1,0,3,0,3,2], [4,3,2,1,0,1],
[1,2,3,2,3,4], [4,3,0,3,0,5] starting at n = 6).

Proof ingredients checked here.
1. Closed form: with W(x) = sum_{m>=1} (-m)^(m-1) x^m/m! (W e^W = x), exp(y W(x)) = sum_m y (y-m)^(m-1) x^m/m!, so
   a(n) = sum_{k=0..n} C(n,k) q^(k(k+1)) (q^k - (n-k))^(n-k-1)   (the k = n term is q^(n^2)).
   Checked against an exact series expansion of the e.g.f. (sympy) for n <= N1.
2. Modulo 2 and 3: if p | q only k = 0 survives, (1-n)^(n-1) mod p, period lcm(p, p-1); if p does not divide q the sum
   is a binomial transform of a function of (n mod 6, k mod 6) plus boundary terms, and (1+x)^n mod (x^6 - 1) over
   F_3 (resp. (x^2 - 1) over F_2) is periodic from n = 3 (resp. 2). Checked: (1+x)^9 = (1+x)^3 mod (x^6-1) over F_3.
3. The patterns for 6 <= n <= N2 from the closed form, and the pattern from n = 3 on.
"""
import sys
from math import comb, factorial
from fractions import Fraction
N1 = int(sys.argv[1]) if len(sys.argv) > 1 else 12
N2 = int(sys.argv[2]) if len(sys.argv) > 2 else 400
def a(n, q):
    s = Fraction(0)
    for k in range(n + 1):
        e = n - k - 1; b = q ** k - (n - k)
        s += comb(n, k) * q ** (k * (k + 1)) * (Fraction(b) ** e if not (b == 0 and e < 0) else 0)
    assert s.denominator == 1
    return int(s)
from sympy import symbols, Rational, factorial as sfac, expand, series
x = symbols('x')
W = sum(Rational((-m) ** (m - 1), sfac(m)) * x ** m for m in range(1, N1 + 1))
ok_series = True
for q in (2, 3, 4, 5):
    A = 0
    for n in range(N1 + 1):
        A += expand((q ** n * x + W) ** n) / sfac(n)
    A = expand(A)
    for n in range(N1 + 1):
        if A.coeff(x, n) * sfac(n) != a(n, q): ok_series = False
DATA = {2: [1, 3, 23, 708, 82677, 39043808, 75384175459], 3: [1, 4, 98, 21901, 45203076, 864855654349],
        4: [1, 5, 287, 274532, 4362420261], 5: [1, 6, 674, 2000229, 153566609748]}
ok_data = all([a(n, q) for n in range(len(DATA[q]))] == DATA[q] for q in DATA)
PAT = {2: [1, 0, 3, 0, 3, 2], 3: [4, 3, 2, 1, 0, 1], 4: [1, 2, 3, 2, 3, 4], 5: [4, 3, 0, 3, 0, 5]}
ok_pat = all(a(n, q) % 6 == PAT[q][(n - 6) % 6] for q in PAT for n in range(6, N2 + 1))
ok_from3 = all(a(n, q) % 6 == PAT[q][(n - 6) % 6] for q in PAT for n in range(3, 6))
def polmul(a_, b_, p, m):
    r = [0] * m
    for i, u in enumerate(a_):
        for j, v in enumerate(b_):
            r[(i + j) % m] = (r[(i + j) % m] + u * v) % p
    return r
def onepx_pow(n, p, m):
    r = [1] + [0] * (m - 1); b = [1, 1] + [0] * (m - 2)
    for _ in range(n): r = polmul(r, b, p, m)
    return r
ok_ring3 = onepx_pow(9, 3, 6) == onepx_pow(3, 3, 6)
ok_ring2 = onepx_pow(2, 2, 2) == [0, 0]
print(f"closed form = series expansion, n<={N1}: {ok_series}; data: {ok_data}")
print(f"patterns mod 6 for 6<=n<={N2}: {ok_pat}; already from n = 3: {ok_from3}; (1+x)^9 = (1+x)^3 mod (x^6-1, 3): {ok_ring3}; (1+x)^2 = 0 mod (x^2-1, 2): {ok_ring2}")
good = ok_series and ok_data and ok_pat and ok_from3 and ok_ring3 and ok_ring2
print("RESULT", "PASS" if good else "FAIL")
sys.exit(0 if good else 1)
