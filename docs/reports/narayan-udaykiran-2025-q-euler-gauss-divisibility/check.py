#!/usr/bin/env python3
"""Checks for S. Narayan and N. Uday Kiran, "On Some Generalisations of Gauss Sequences", arXiv:2511.19503v5,
Conjecture 4.10 (Conjecture 4.9 in v1-v3).

Definition 4.1: (a_n(q)) in Z[q] is q-Euler-Gauss if prod_{d|n, mu(d)=1} a_{n/d}(q^d) = prod_{d|n, mu(d)=-1} a_{n/d}(q^d)
mod [n]_q for all n >= 1. Conjecture 4.10: then [n]_q divides S_n(q) A_n(q), S_n = sum_{d|n} mu(d) a_{n/d}(q^d),
A_1 = 1, A_n = prod_{d|n, d>1, mu(d)=1} a_{n/d}(q^d) * prod_{d|n, mu(d)=-1} A_{n/d}(q^d).
Part 1 (exact, Z[q]): Examples 4.11 and 4.14 satisfy (4.1) and fail the q-Gauss congruence exactly where the source says
(fidelity); Example 4.13 fails (4.1) (negative control); the criterion holds for both examples for n <= 60.
Part 2 (finite field F_p, p = 2 lcm(1..60) + 1): random sequences of values at all roots of unity of order <= 60 that
satisfy (4.1) at every nontrivial root; A_n S_n vanishes at every nontrivial n-th root, n <= 60, while S_n alone and
B_n S_n do not (power of the test). Exits 1 on failure.
"""
import math, random, sys
import sympy as sp
from sympy import mobius, divisors, factorint, isprime

ok = True
q = sp.symbols('q')
N = 60
def qint(n): return sp.Poly(sum(q**j for j in range(n)), q, domain='ZZ')
def ex414(n):
    if n == 1: return sp.Poly(q + 1, q, domain='ZZ')
    if n % 2: return sp.Poly(2, q, domain='ZZ')
    m = n // 2
    if len(factorint(m)) <= 1: return sp.Poly(1 - q**m, q, domain="ZZ")
    Om = sum(factorint(n).values())
    return sp.Poly((2**(Om - 3) - 1) * q**m + (2**(Om - 3) + 1), q, domain='ZZ')
def ex411(n):
    f = factorint(n)
    if n % 12: return sp.Poly(q**n - 1, q, domain='ZZ')
    k = f.get(2, 0)
    return sp.Poly(q**(n // (2**k * 3)), q, domain='ZZ')
def sub(P, d): return sp.Poly(P.as_expr().subs(q, q**d), q, domain='ZZ')
def ex413(n):
    if n % 2: return sp.Poly(1, q, domain='ZZ')
    k = factorint(n).get(2, 0)
    return sp.Poly(sum(q**j for j in range(2)), q, domain='ZZ').compose(sp.Poly(q**(2**(k - 1)), q, domain='ZZ')) + sp.Poly(1, q, domain='ZZ')
def run(a, name):
    seq = {n: a(n) for n in range(1, N + 1)}
    A = {1: sp.Poly(1, q, domain='ZZ')}
    eg_fail = []; conj_fail = []; gauss_fail = []
    for n in range(1, N + 1):
        Qn = qint(n)
        Pp = sp.Poly(1, q, domain='ZZ'); Pm = sp.Poly(1, q, domain='ZZ'); S = sp.Poly(0, q, domain='ZZ')
        for d in divisors(n):
            m = int(mobius(d))
            t = sub(seq[n // d], d)
            if m == 1: Pp = (Pp * t).rem(Qn)
            if m == -1: Pm = (Pm * t).rem(Qn)
            S = S + m * t
        if not (Pp - Pm).rem(Qn).is_zero: eg_fail.append(n)
        if not S.rem(Qn).is_zero: gauss_fail.append(n)
        if n > 1:
            B = sp.Poly(1, q, domain='ZZ')
            for d in divisors(n):
                if d > 1 and mobius(d) == 1: B = B * sub(seq[n // d], d)
            M = sp.Poly(1, q, domain='ZZ')
            for d in divisors(n):
                if mobius(d) == -1: M = M * sub(A[n // d], d)
            A[n] = B * M
            if not (A[n] * S).rem(Qn).is_zero: conj_fail.append(n)
    print(f"{name}: (4.1) fails at {eg_fail}; q-Gauss fails at {gauss_fail}; criterion [n]_q | A_n S_n fails at {conj_fail} (n <= {N}); deg A_N = {A[N].degree()}")
    return eg_fail, gauss_fail, conj_fail

e414 = run(ex414, 'Example 4.14')
want414 = [n for n in range(2, 61, 2) if len(factorint(n // 2)) > 1]
f411 = run(ex411, 'Example 4.11')
e413 = run(ex413, 'Example 4.13 (negative control)')
print(f"fidelity: Example 4.14 fails q-Gauss exactly at the even n <= 60 not of the form 2p^k: {e414[1] == want414}; Example 4.11 exactly at 2^k 3: {f411[1] == [12, 24, 48]}; Example 4.13 is not q-Euler-Gauss: {len(e413[0]) > 0}")
ok &= e414[0] == [] and e414[2] == [] and e414[1] == want414 and f411[0] == [] and f411[2] == [] and f411[1] == [12, 24, 48] and len(e413[0]) > 0
N = 60
L = 1
for k in range(1, N + 1): L = L * k // math.gcd(L, k)
k = 1
while not isprime(k * L + 1): k += 1
p = k * L + 1
fac = set(factorint(k)) | set(factorint(L))
g = 2
while any(pow(g, (p - 1) // q, p) == 1 for q in fac): g += 1
G = pow(g, (p - 1) // L, p)
order = lambda e: L // math.gcd(e, L)
# enumerate exponents with order <= N: e = (L/r)*j, gcd(j, r) = 1
E = sorted({(L // r) * j % L for r in range(1, N + 1) for j in range(r) if math.gcd(j, r) == 1})
mu = {n: int(mobius(n)) for n in range(1, N + 1)}
divs = {n: divisors(n) for n in range(1, N + 1)}
def trial(rng, pzero):
    v = {}
    for n in range(1, N + 1):
        v[n] = {}
        for e in E:
            o = order(e)
            if o > 1 and n % o == 0:
                B = 1; Pm = 1
                for d in divs[n]:
                    if d > 1 and mu[d] == 1: B = B * v[n // d][e * d % L] % p
                    if mu[d] == -1: Pm = Pm * v[n // d][e * d % L] % p
                if B: v[n][e] = Pm * pow(B, p - 2, p) % p
                elif Pm: return None
                else: v[n][e] = 0 if rng.random() < pzero else rng.randrange(1, p)
            else:
                v[n][e] = 0 if rng.random() < pzero else rng.randrange(1, p)
    return v
def conclusions(v):
    A = {1: {e: 1 for e in E}}
    fails = {'AS': 0, 'S': 0, 'BS': 0}; checked = 0
    for n in range(2, N + 1):
        A[n] = {}
        for e in E:
            B = 1
            for d in divs[n]:
                if d > 1 and mu[d] == 1: B = B * v[n // d][e * d % L] % p
            M = 1
            for d in divs[n]:
                if mu[d] == -1: M = M * A[n // d][e * d % L] % p
            A[n][e] = B * M % p
        for e in E:
            o = order(e)
            if o > 1 and n % o == 0:
                S = sum(mu[d] * v[n // d][e * d % L] for d in divs[n]) % p
                B = 1
                for d in divs[n]:
                    if d > 1 and mu[d] == 1: B = B * v[n // d][e * d % L] % p
                checked += 1
                fails['AS'] += (A[n][e] * S % p != 0)
                fails['S'] += (S != 0)
                fails['BS'] += (B * S % p != 0)
    return fails, checked
rng = random.Random(20261023)
tot = {'AS': 0, 'S': 0, 'BS': 0}; ok_trials = 0; chk = 0; restarts = 0
for t in range(60):
    pz = [0.2, 0.4, 0.6][t % 3]
    v = None
    while v is None:
        v = trial(rng, pz); restarts += v is None
    f, c = conclusions(v); ok_trials += 1; chk += c
    for kk in tot: tot[kk] += f[kk]
print(f"p = {k}*lcm(1..{N})+1; trials {ok_trials} (restarts {restarts}); root checks {chk}; failures of A_n S_n = 0: {tot['AS']}; of S_n = 0: {tot['S']}; of B_n S_n = 0: {tot['BS']}")
ok &= tot['AS'] == 0 and tot['S'] > 0 and tot['BS'] > 0
print("ALL_CHECKS_OK" if ok else "CHECK_FAIL")
sys.exit(0 if ok else 1)
