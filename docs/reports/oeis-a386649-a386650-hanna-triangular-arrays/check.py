"""OEIS A386649 / A386650 (P. D. Hanna): products of central trinomial (resp. quadrinomial)
coefficients versus R. H. Hardin's counts A214589 / A214590 of n X n X n triangular arrays
with entries 0..d (d = 2, resp. 3) whose rows (lengths 1..n) all have the same average.

Conjectures: A386649(n) = A214589(n) - 2 and 2*A386650(n) = A214590(n) - 2 for n >= 1.
The script (1) counts arrays by brute force for small n, (2) counts them for n <= 13 by
summing, over the admissible common averages, the product of row counts, (3) compares
with both published data and with the product formulas.
"""
import itertools, sys
from fractions import Fraction

A214589 = [3,5,23,401,20351,2869211,1127599139,1248252244661,3918263795984603,
           35080215765450132755,899912775031092255512711,66403663756769266442027284403,
           14140062564030204365431731967633343]
A214590 = [4,8,62,1862,187862,63120962,71200442882,272199293130242,3552472974642762242,
           159161446682919676638722,24593308418551382594861736962,13156534644821921838477614251069442]

def rows_with_sum(length, d, s):
    # number of sequences in {0..d}^length with sum s
    c = [1] + [0] * (d * length)
    for _ in range(length):
        nc = [0] * (d * length + 1)
        for t, v in enumerate(c):
            if v:
                for e in range(d + 1):
                    if t + e < len(nc): nc[t + e] += v
        c = nc
    return c[s] if 0 <= s < len(c) else 0

def coeff(d, n):  # [x^n] (1 + x + ... + x^d)^n
    return rows_with_sum(n, d, n)

def count_formula(n, d):
    tot = 0
    for r in range(d + 1):  # the first row fixes the average to an integer 0..d
        prod = 1
        for i in range(1, n + 1):
            prod *= rows_with_sum(i, d, r * i)
        tot += prod
    return tot

def count_brute(n, d):
    tot = 0
    rows = [list(itertools.product(range(d + 1), repeat=i)) for i in range(1, n + 1)]
    for arr in itertools.product(*rows):
        avgs = {Fraction(sum(r), len(r)) for r in arr}
        tot += len(avgs) == 1
    return tot

ok = True
for d, data in ((2, A214589), (3, A214590)):
    for n in range(1, 4 if d == 3 else 5):
        if count_brute(n, d) != data[n - 1]: ok = False; print("brute mismatch", d, n)
    for n in range(1, len(data) + 1):
        f = count_formula(n, d)
        prod = 1
        for k in range(n + 1): prod *= coeff(d, k)
        want = prod + 2 if d == 2 else 2 * prod + 2
        if not (f == data[n - 1] == want): ok = False; print("mismatch", d, n, f, data[n - 1], want)
print(f"brute force (d=2: n<=4, d=3: n<=3), average-sum count and product formulas agree with the data: {ok}")
print("RESULT", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
