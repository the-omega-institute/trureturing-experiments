"""OEIS A384531 (W. Schulte) and A367865 (W. I. Hurt): two multiplicative gcd-sum identities.

A384531: a(n) multiplicative with a(p^e) = ((2e+1)p - 2e) p^(e-1).
  Conjecture: a(n) = Sum_{i=1..n} gcd(i, n) * gcd(i+1, n).
A367865: a(n) = Sum_{d|n} d * phi(d) * mu(d)^2.
  Conjecture: a(n) = (rad(n)/n) * Sum_{k=1..n} rad(n*k)/rad(k).
For n <= N the script checks both conjectures exactly (rational arithmetic for A367865), the
published data, and the proof's local facts: at most one of gcd(i, p^e), gcd(i+1, p^e) exceeds 1,
and rad(n*k)/rad(k) = rad(n)/rad(gcd(n, k)).
"""
import sys
from math import gcd
from fractions import Fraction
from sympy import factorint, totient, mobius

N = int(sys.argv[1]) if len(sys.argv) > 1 else 1500
D1 = [1,4,7,12,13,28,19,32,33,52,31,84,37,76,91,80,49,132,55,156,133,124]
D2 = [1,3,7,3,21,21,43,3,7,63,111,21,157,129,147,3,273,21,343,63,301,333]

def rad(n):
    r = 1
    for p in factorint(n): r *= p
    return r

def a1(n):
    r = 1
    for p, e in factorint(n).items(): r *= ((2*e+1)*p - 2*e) * p**(e-1)
    return r

def a2(n):
    return sum(d * totient(d) * mobius(d)**2 for d in range(1, n + 1) if n % d == 0)

ok1 = all(a1(n) == sum(gcd(i, n) * gcd(i + 1, n) for i in range(1, n + 1)) for n in range(1, N + 1))
ok2 = all(Fraction(a2(n)) == Fraction(rad(n), n) * sum(Fraction(rad(n * k), rad(k)) for k in range(1, n + 1)) for n in range(1, N // 3 + 1))
okd = [a1(n) for n in range(1, len(D1) + 1)] == D1 and [a2(n) for n in range(1, len(D2) + 1)] == D2
okl = all(not (gcd(i, q) > 1 and gcd(i + 1, q) > 1) for q in (2, 4, 8, 9, 27, 25, 49, 121) for i in range(q)) and \
      all(Fraction(rad(n * k), rad(k)) == Fraction(rad(n), rad(gcd(n, k))) for n in range(1, 120) for k in range(1, 120))
print(f"A384531 conjecture for n <= {N}: {ok1}; A367865 conjecture for n <= {N // 3}: {ok2}; data: {okd}; local facts: {okl}")
good = ok1 and ok2 and okd and okl
print("RESULT", "PASS" if good else "FAIL")
sys.exit(0 if good else 1)
