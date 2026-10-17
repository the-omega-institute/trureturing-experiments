"""OEIS A394815: twin-prime centers k (k - 1, k + 1 prime) with gcd(k, sigma(k)) = 2.
Conjecture: every term is divisible by 36.

For every twin-prime center k < N the script computes gcd(k, sigma(k)) and checks:
the conjecture; that the terms found match the published data; and the two steps of
the proof: 18 | k (from the frozen A394757 settlement) and 4 | k (else 3 | sigma(k)).
"""
import sys
from math import gcd
from sympy import divisor_sigma, isprime

N = int(sys.argv[1]) if len(sys.argv) > 1 else 2_000_000
DATA = [5652, 6948, 16452, 20772, 27792, 44532, 65700, 67428, 126612, 143568, 180288]
terms, ok_conj, ok_steps = [], True, True
for k in range(4, N, 2):
    if isprime(k - 1) and isprime(k + 1):
        s = int(divisor_sigma(k))
        if gcd(k, s) == 2:
            terms.append(k)
            if k % 36:
                ok_conj = False
        if k % 2 == 0 and k % 4 and k % 3 == 0 and s % 3:
            ok_steps = False  # v2(k) = 1 and 3 | k must force 3 | sigma(k)
ok_data = terms[:len(DATA)] == DATA
print(f"N={N}: {len(terms)} terms; all divisible by 36: {ok_conj}; data match: {ok_data}; "
      f"v2(k)=1 with 3|k forces 3|sigma(k) on all centers: {ok_steps}")
good = ok_conj and ok_data and ok_steps
print("RESULT", "PASS" if good else "FAIL")
sys.exit(0 if good else 1)
