"""OEIS A372761 (x = n+4), A372763 (x = n+5), A363347 (x = -4), A363102 (x = -2): denominators of
K_n(x) = 1/(2-3/(3-4/(4-...(n-1)-n/x))), n >= 3.

Checked statements:
- A372761, Conjecture 1: except for 4, every denominator is 1 or a prime.
- A372763, Conjecture 2: except for 2 and 3, every prime is a denominator at exactly one n. (Refuted: 11 occurs at n = 10 and 21.)
- A363347, Conjecture 3: every prime denominator other than 5 occurs at exactly two n.
- A363102, Conjecture 2: every prime denominator occurs at exactly two or three n.

The script (1) evaluates the nested fraction exactly for 3 <= n <= N1, checks that no intermediate
quotient is zero, compares with the published data, and checks the closed form
K_n(x) = (a(n-3)x - n a(n-4))/((n-1)x - n(n-2)) (a = A051403) and den = |D|/gcd(|D|, n!/2) with
D = (n-1)x - n(n-2); (2) computes den = prod_p p^max(0, v_p(D) - v_p(n!/2)) for 3 <= n <= N2;
(3) checks the four statements and the proofs' predictions for every prime p <= N2/2, all of whose
occurrences lie below 2p <= N2:
- A372763: 2 and 3 never occur, 11 occurs at n = 10 and 21, every other prime 5 <= p <= N2/2 exactly once;
- A363347: 5 occurs only at n = 4, every other occurring prime exactly twice, at the two roots
  n < p of (n+1)^2 = 5 (mod p);
- A363102: every occurring prime occurs at the two roots n < p of n^2 = 2 (mod p), and a third
  time exactly when 2p^2 + 2 is a square (at n = sqrt(2p^2 + 2)).
"""
import sys
from fractions import Fraction
from math import gcd, factorial, isqrt
from collections import defaultdict
from sympy import isprime, factorint, primerange

N1 = int(sys.argv[1]) if len(sys.argv) > 1 else 300
N2 = int(sys.argv[2]) if len(sys.argv) > 2 else 20000
CASES = {'A372761': (lambda n: n + 4, [11,4,7,13,31,1,41,23,17,1,61,1,71,19,1,43,1,1,101,53,37,29,1,1,131,1]),
         'A372763': (lambda n: n + 5, [13,19,5,31,37,43,7,11,61,67,73,79,17,1,97,103,109,23,11,127,1,139,29]),
         'A363347': (lambda n: -4, [11,5,31,11,59,19,19,29,139,41,191,1,251,71,29,89,79,109,479,131,571]),
         'A363102': (lambda n: -2, [7,7,23,17,47,31,79,7,17,71,167,97,223,127,41,23,359,199,439,241,31])}
a = {-1: 0, 0: 1}
for m in range(1, N1 + 2): a[m] = (m + 2) * (a[m - 1] - a[m - 2])

def nested(n, x):
    t = Fraction(x)
    for k in range(n - 1, 1, -1):
        if t == 0: return None
        t = k - Fraction(k + 1) / t
    if t == 0: return None
    return 1 / t

def vfact(n, p):
    s, q = 0, p
    while q <= n: s += n // q; q *= p
    return s

def den(n, x):
    D = abs((n - 1) * x - n * (n - 2)); r = 1
    for p, e in factorint(D).items():
        f = vfact(n, p) - (1 if p == 2 else 0)
        if e > f: r *= p ** (e - f)
    return r

good = True
dens = {}
for name, (xf, data) in CASES.items():
    exact, ok_closed, ok_nonzero = [], True, True
    for n in range(3, N1 + 1):
        x = xf(n); v = nested(n, x)
        if v is None: ok_nonzero = False; exact.append(None); continue
        D = (n - 1) * x - n * (n - 2)
        if v != Fraction(a[n - 3] * x - n * a[n - 4], D): ok_closed = False
        if v.denominator != abs(D) // gcd(abs(D), factorial(n) // 2): ok_closed = False
        exact.append(v.denominator)
    ok_data = exact[:len(data)] == data
    d = {n: den(n, xf(n)) for n in range(3, N2 + 1)}
    ok_val = all(d[n] == exact[n - 3] for n in range(3, N1 + 1))
    dens[name] = d
    good = good and ok_nonzero and ok_closed and ok_data and ok_val
    print(f"{name}: nested nonzero n<={N1}: {ok_nonzero}; closed form + den formula: {ok_closed}; "
          f"data: {ok_data}; valuation formula = exact n<={N1}: {ok_val}")

def occ(d):
    o = defaultdict(list)
    for n, v in d.items():
        if v != 1: o[v].append(n)
    return o

P = list(primerange(2, N2 // 2 + 1))
# A372761, Conjecture 1
bad = [(n, v) for n, v in dens['A372761'].items() if v != 1 and not isprime(v)]
ok = bad == [(4, 4)]
good = good and ok
print(f"A372761 Conjecture 1: non-prime denominators > 1 for n<={N2}: {bad}; holds: {ok}")
# A372763, Conjecture 2
o = occ(dens['A372763'])
pred = {p: ([] if p in (2, 3) else [10, 21] if p == 11 else None) for p in P}
dev = [(p, o.get(p, [])) for p in P if (pred[p] is not None and o.get(p, []) != pred[p]) or (pred[p] is None and len(o.get(p, [])) != 1)]
ok = dev == []
good = good and ok
print(f"A372763 Conjecture 2: occurrences of 11: {o.get(11)}; occurrences of 5, 7: {o.get(5)}, {o.get(7)}; "
      f"deviations from the corrected classification for p<={N2//2}: {dev}; refuted at p = 11: {o.get(11) == [10, 21]}")
# A363347, Conjecture 3
o = occ(dens['A363347'])
dev = []
for p in P:
    ns = o.get(p, [])
    if p == 5: want = [4]
    else:
        roots = [n for n in range(3, p) if ((n + 1) ** 2 - 5) % p == 0]
        want = roots if len(roots) == 2 else []
    if ns != want: dev.append((p, ns, want))
cnt = defaultdict(int)
for p in P: cnt[len(o.get(p, []))] += 1
ok = dev == []
good = good and ok
print(f"A363347 Conjecture 3: occurrence counts for p<={N2//2}: {dict(cnt)}; 5 at {o.get(5)}; deviations from the roots prediction: {dev[:5]}; holds: {ok}")
# A363102, Conjecture 2
o = occ(dens['A363102'])
dev, triples = [], []
for p in P:
    ns = o.get(p, [])
    roots = [n for n in range(3, p) if (n * n - 2) % p == 0]
    want = roots if len(roots) == 2 else []
    s = 2 * p * p + 2; r = isqrt(s)
    if want and r * r == s: want = want + [r]; triples.append((p, r))
    if ns != want: dev.append((p, ns, want))
cnt = defaultdict(int)
for p in P: cnt[len(o.get(p, []))] += 1
ok = dev == [] and all(c in (0, 2, 3) for c in cnt)
good = good and ok
print(f"A363102 Conjecture 2: occurrence counts for p<={N2//2}: {dict(cnt)}; three times (p, third n): {triples}; deviations: {dev[:5]}; holds: {ok}")
print("RESULT", "PASS" if good else "FAIL")
sys.exit(0 if good else 1)
