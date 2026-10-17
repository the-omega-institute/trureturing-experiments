"""OEIS A400454 / A400825 (G. Wright): for N = (2^k - 1)*2^(k + j) + 1 with
j = 3 (A400454) or j = 4 (A400825), the conjecture "if 2^((N-1)/4) == -1 (mod N)
then N is prime".

For 1 <= k <= K the script records whether the hypothesis holds and whether N is
prime (sympy.isprime), checks the implication, checks the size bound
N < (2^(k+j-1) + 1)^2 used in the proof, and compares the primes found with the
published terms.
"""
import sys
from sympy import isprime

K = int(sys.argv[1]) if len(sys.argv) > 1 else 400
DATA = {3: [1,2,3,5,6,10,11,14,16,18,22,26,34,62,107,135,136,270,271,318,371],
        4: [2,6,10,12,16,30,60,64,94,118,120,138,178,200,226]}
good = True
for j in (3, 4):
    hyp, primes, viol, bound_ok = [], [], [], True
    for k in range(1, K + 1):
        N = (2 ** k - 1) * 2 ** (k + j) + 1
        h = pow(2, (N - 1) // 4, N) == N - 1
        p = isprime(N)
        if h: hyp.append(k)
        if p: primes.append(k)
        if h and not p: viol.append(k)
        if not N < (2 ** (k + j - 1) + 1) ** 2: bound_ok = False
    pref = [k for k in primes if k <= max(DATA[j])]
    ok = not viol and bound_ok and pref == DATA[j]
    good = good and ok
    print(f"j={j}, 1<=k<={K}: hypothesis holds for k in {hyp[:12]}... ({len(hyp)} values); "
          f"violations {viol}; size bound {bound_ok}; primes match data {pref == DATA[j]}; "
          f"hypothesis-true primes / all primes = {len(hyp)}/{len(primes)}")
print("RESULT", "PASS" if good else "FAIL")
sys.exit(0 if good else 1)
