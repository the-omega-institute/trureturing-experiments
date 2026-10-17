"""OEIS A350532 (P. Kagey): T(n,k) = number of degree-n polynomials over GF(2) of the form
f(x)^m (m > 1) with exactly k nonzero terms. Conjectures:
  (1) T(n,n+1) = 1 iff n is a Mersenne prime;  (2) T(2n,2) = n;  (3) T(2n,3) = (n^2-n)/2.

A polynomial over GF(2) is f^m for some m > 1 iff the gcd of the multiplicities in its
factorization is > 1 (unique factorization). The script tests this with sympy's GF(2)
factorization and checks: the published rows n <= 10 (all polynomials); T(n, n+1) for n <= 70
(only the all-ones polynomial has n+1 terms); T(2n,2) for n <= 40 and T(n,3) for n <= 40
(all binomials and trinomials); and the explicit counterexample to (1):
1 + x + ... + x^15 = (1 + x)^15 = ((1 + x)^5)^3 with 15 = 3 * 5 not prime.
"""
import sys
from math import gcd
from itertools import combinations
from sympy import Poly, symbols, isprime
x = symbols('x')
DATA = [1,0,0,1,1,0,1,0,0,1,1,2,1,0,0,1,0,0,1,0,0,1,3,3,2,1,0,0,1,0,0,0,0,0,0,1,1,4,6,4,1,0,0,0,0,1,0,0,4,1,0,2,0,0,0,1,5,10,11,5,1,0,0,1,0,0]

def is_power(exps):
    if not exps: return False
    p = Poly(sum(x**e for e in exps), x, modulus=2)
    if p.degree() < 1: return False
    _, facs = p.factor_list()
    g = 0
    for _, e in facs: g = gcd(g, e)
    return g > 1

def row(n):
    r = [0] * (n + 2)
    for mask in range(1 << n):           # coefficients of x^0..x^(n-1); x^n always present
        exps = [i for i in range(n) if mask >> i & 1] + [n]
        if is_power(exps): r[len(exps)] += 1
    return r[1:]
tri = []
for n in range(0, 11):
    tri += [1] if n == 0 else row(n)
ok_data = tri == DATA
c1 = [n for n in range(1, 71) if is_power(list(range(n + 1)))]
mers = [n for n in range(1, 71) if isprime(n) and (n + 1) & n == 0]
pow2 = [n for n in range(3, 71) if (n + 1) & n == 0]
ok_c1_refuted = c1 == pow2 and 15 in c1 and 15 not in mers
ok_c2 = all(sum(is_power([b, 2 * n]) for b in range(2 * n)) == n for n in range(1, 41))
def t3(n): return sum(is_power([b, a, n]) for a, b in combinations(range(n), 2))
ok_c3 = all(t3(2 * n) == (n * n - n) // 2 for n in range(1, 21)) and all(t3(n) == 0 for n in range(1, 40, 2))
cex = Poly(sum(x**i for i in range(16)), x, modulus=2) == Poly((1 + x)**15, x, modulus=2) == Poly(((1 + x)**5)**3, x, modulus=2)
print(f"rows n<=10 match data: {ok_data}")
print(f"n<=70 with T(n,n+1)=1: {c1}; Mersenne primes <=70: {mers}; equals {{2^j-1 : j>=2}}: {c1 == pow2}")
print(f"counterexample 1+x+...+x^15 = (1+x)^15 = ((1+x)^5)^3 over GF(2): {cex}; 15 prime: {isprime(15)}")
print(f"T(2n,2) = n for n<=40: {ok_c2}; T(2n,3) = (n^2-n)/2 for n<=20 and T(odd,3) = 0 below 40: {ok_c3}")
good = ok_data and ok_c1_refuted and cex and ok_c2 and ok_c3
print("RESULT", "PASS" if good else "FAIL")
sys.exit(0 if good else 1)
