"""OEIS A356360, A370726, A372763 (x = n+1, n+3, n+5) and A362086, A363482 (x = -3, -5):
denominators of K_n(x) = 1/(2-3/(3-4/(4-...(n-1)-n/x))). Conjectures: the denominators are
1 or primes (A356360, A370726, A372763), resp. 1 or primes except the single value 9 (A362086)
and 49 (A363482).

The script (1) evaluates the nested fraction exactly for 3 <= n <= N1, checks that no
intermediate quotient is zero, and compares with the published data; (2) checks the proof's
closed form K_n(x) = (a(n-3)x - n a(n-4))/((n-1)x - n(n-2)) (a = A051403) and the denominator
formula den = |D|/gcd(|D|, n!/2) with D = (n-1)x - n(n-2); (3) checks the conjectures with the
formula for 3 <= n <= N2, recording every composite denominator.
"""
import sys
from fractions import Fraction
from math import gcd, factorial
from sympy import isprime

N1 = int(sys.argv[1]) if len(sys.argv) > 1 else 300
N2 = int(sys.argv[2]) if len(sys.argv) > 2 else 4000
CASES = {'A356360': (lambda n: n + 1, [5,7,3,11,13,1,17,19,1,23,1,1,29,31,1,1,37,1,41,43,1,47,1,1,53,1,1,59], set()),
         'A370726': (lambda n: n + 3, [3,13,17,7,5,29,11,37,41,1,7,53,19,61,1,23,73,1,1,1,89,31,97,101,1], set()),
         'A372763': (lambda n: n + 5, [13,19,5,31,37,43,7,11,61,67,73,79,17,1,97,103,109,23,11,127,1,139,29], set()),
         'A362086': (lambda n: -3, [3,17,9,13,53,23,29,107,43,17,179,23,79,269,101,113,29,139,1,503,61], {9}),
         'A363482': (lambda n: -5, [13,23,7,49,13,83,103,5,149,1,29,233,53,23,67,373,59,1,499,109,593], {49})}
a = {-1: 0, 0: 1}
for m in range(1, N2 + 2): a[m] = (m + 2) * (a[m - 1] - a[m - 2])

def nested(n, x):
    t = Fraction(x)
    for k in range(n - 1, 1, -1):
        if t == 0: return None
        t = k - Fraction(k + 1) / t
    if t == 0: return None
    return 1 / t

good = True
for name, (xf, data, allowed) in CASES.items():
    dens, ok_closed, ok_nonzero = [], True, True
    for n in range(3, N1 + 1):
        x = xf(n); v = nested(n, x)
        if v is None: ok_nonzero = False; dens.append(None); continue
        D = (n - 1) * x - n * (n - 2)
        if v != Fraction(a[n - 3] * x - n * a[n - 4], D): ok_closed = False
        if v.denominator != abs(D) // gcd(abs(D), factorial(n) // 2): ok_closed = False
        dens.append(v.denominator)
    ok_data = dens[:len(data)] == data
    bad = []
    for n in range(3, N2 + 1):
        x = xf(n); D = abs((n - 1) * x - n * (n - 2))
        d = D // gcd(D, factorial(n) // 2)
        if d != 1 and not isprime(d): bad.append((n, d))
    ok_conj = all(d in allowed for _, d in bad) and len(bad) == len(allowed)
    good = good and ok_data and ok_closed and ok_nonzero and ok_conj
    print(f"{name}: nested nonzero n<={N1}: {ok_nonzero}; closed form + den formula: {ok_closed}; data: {ok_data}; "
          f"composite denominators n<={N2}: {bad}; conjecture: {ok_conj}")
print("RESULT", "PASS" if good else "FAIL")
sys.exit(0 if good else 1)
