"""McKean, Toda primes (arXiv:2511.19744v1) Conjecture 5.6: counterexamples.

T(n) = {odd primes p : p-1 | 4n, gcd(p, 4n/(p-1)) = 1}, t(n) = |T(n)|; D_k = denominator of B_k (direct for k <= 5000, von Staudt-Clausen above, cross-checked for k <= 400);
F(d) = min{k even > 0 : D_k = d}.
(i)  F(D_{4a}) = 4a and D_{4am} = D_{4a}  =>  T(am) = T(a)
(ii) D_{4am} != D_{4a}                      =>  t(am) >= t(a)
Usage: PYTHONPATH=<sympy> python3 mckean-c56.py   (exit 0 iff all readings hold)
"""
import sys
from math import gcd
from sympy import bernoulli, isprime, primerange
def T(n):
    return {p for p in primerange(3, 4*n+2) if (4*n) % (p-1) == 0 and gcd(p, 4*n//(p-1)) == 1}
from sympy import divisors
def D_vsc(k):  # von Staudt-Clausen: product of primes p with p-1 | k
    r = 1
    for e in divisors(k):
        if isprime(e + 1): r *= e + 1
    return r
def D(k):
    return D_vsc(k) if k > 5000 else bernoulli(k).q
def F(d, bound):
    for k in range(2, bound+1, 2):
        if D_vsc(k) == d: return k
    return None
ok = True
# cross-check von Staudt-Clausen against direct Bernoulli denominators
for k in range(2, 401, 2):
    assert bernoulli(k).q == D_vsc(k), k
print("VSC cross-check k<=400: ok")
# (ii) witness a = 11, m = 5
a, m = 11, 5
r = (D(4*a), D(4*a*m), F(D(4*a), 4*a), T(a), T(a*m))
print("(ii) a=11 m=5: D44, D220, F(D44), T(11), T(55) =", r)
ok &= r[0] == 690 and r[1] == 7590 and r[2] == 44 and len(r[4]) < len(r[3])
# (i) witness a = 23, m = 47
a, m = 23, 47
r = (D(4*a), D(4*a*m), F(D(4*a), 4*a), T(a), T(a*m))
print("(i) a=23 m=47: D92, D4324, F(D92), T(23), T(1081) =", r)
ok &= r[0] == 1410 and r[1] == 1410 and r[2] == 92 and r[3] != r[4]
# family: a prime, a = 23 mod 30, p = 2a+1 prime
fam = [a for a in primerange(23, 2000) if a % 30 == 23 and isprime(2*a+1)]
for a in fam[:8]:
    p = 2*a + 1
    good = D(4*a) == 30*p and D(4*a*p) == 30*p and F(30*p, 4*a) == 4*a and T(a) == {3, 5, p} and T(a*p) == {3, 5}
    print(f"family a={a} p={p}: {good}")
    ok &= good
print("family members below 2000:", fam)
sys.exit(0 if ok else 1)
