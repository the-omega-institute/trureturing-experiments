"""Checks for OEIS A368959 (Kotesovec): Conjecture 1 (A368959 is a subset of the powerful numbers A001694) and
Conjecture 2 (A368959 contains A338325 and 243, 972, 1944, 2401), together with the exact counts used by the proof.
For s >= max exponent of m, k^s mod m is periodic in s with period lambda(m), so s in [2, e + lambda(m) + 1] covers all s >= 2."""
import cmath, math
from sympy import factorint
from sympy.functions.combinatorial.numbers import reduced_totient

def counts(m, s):
    c = [0] * m
    for x in range(m):
        c[pow(x, s, m)] += 1
    return c

def S(m, s):
    return sum(cnt * cmath.exp(2j * math.pi * j / m) for j, cnt in enumerate(counts(m, s)) if cnt)

def s_range(m):
    f = factorint(m); e = max(f.values()) if f else 1
    lam = int(reduced_totient(m)) if m > 1 else 1
    return range(2, e + lam + 2)

def in_a368959(m):
    return all(abs(S(m, s)) > 1e-6 for s in s_range(m))

def powerful(m):
    return all(v >= 2 for v in factorint(m).values())

def in_a338325(m):
    return all(v in (2, 3) for v in factorint(m).values())

ok = True
data = [1, 4, 8, 9, 25, 27, 36, 49, 72, 100, 108, 121, 125, 169, 196, 200, 216, 225, 243, 289]
brute = [m for m in range(1, 290) if in_a368959(m)]
ok &= brute == data
print("A368959 data reproduced (first 20 terms):", brute == data)

c2 = all(in_a368959(m) for m in range(1, 3001) if in_a338325(m)) and all(in_a368959(m) for m in (243, 972, 1944, 2401))
ok &= c2
print("Conjecture 2: A338325 terms <= 3000 and 243, 972, 1944, 2401 are in A368959:", c2)

c1 = True
for m in range(2, 3001):
    if not powerful(m):
        p = next(q for q, v in factorint(m).items() if v == 1)
        c1 &= abs(S(m, p)) < 1e-6
print("Conjecture 1: every non-powerful m <= 3000 has S_p(m) = 0 for a prime p dividing m exactly once:", c1)
ok &= c1

def lemma2(p, a, s):
    q = p ** a; h = p ** (a - 1)
    c0 = p ** (a - math.ceil(a / s))
    if (a - 1) % s:
        ch = 0
    else:
        t = (a - 1) // s; ch = math.gcd(s, p - 1) * p ** ((s - 1) * t)
    c = counts(q, s)
    return c[0] == c0 and c[h] == ch
l2 = all(lemma2(p, a, s) for p in (2, 3, 5, 7, 11, 13) for a in range(2, 6) for s in range(2, 13) if p ** a <= 200000)
ok &= l2
print("Lemma 2 exact counts for p <= 13, 2 <= a <= 5, 2 <= s <= 12:", l2)

print("RESULT", "PASS" if ok else "FAIL")
