"""OEIS A394407: k with no smaller m having phi(m)/tau(m) = phi(k)/tau(k).

Checks (1) the first 64 published terms, (2) that every prime p != 3 below N
is a term, and (3) the inequality phi(m)/tau(m) <= (m - 1)/2 for 2 <= m < N
used in the proof.
"""
from fractions import Fraction
import sys

N = int(sys.argv[1]) if len(sys.argv) > 1 else 100000
DATA = [1, 2, 4, 5, 7, 11, 13, 14, 16, 17, 19, 20, 22, 23, 25, 27, 29, 31, 32, 34, 37, 41, 43,
        44, 46, 47, 53, 54, 55, 58, 59, 61, 62, 64, 65, 67, 68, 71, 73, 79, 80, 81, 83, 85, 86,
        89, 92, 94, 97, 100, 101, 103, 106, 107, 109, 112, 113, 115, 116, 118, 119, 121, 125, 127]

phi = list(range(N))
for i in range(2, N):
    if phi[i] == i:
        for j in range(i, N, i):
            phi[j] -= phi[j] // i
tau = [0] * N
for d in range(1, N):
    for j in range(d, N, d):
        tau[j] += 1
is_prime = [False, False] + [phi[i] == i - 1 for i in range(2, N)]

seen, terms = set(), []
for k in range(1, N):
    r = Fraction(phi[k], tau[k])
    if r not in seen:
        terms.append(k)
    seen.add(r)
term_set = set(terms)

ok_data = terms[:len(DATA)] == DATA
missing = [p for p in range(2, N) if is_prime[p] and p != 3 and p not in term_set]
ok_ineq = all(2 * phi[m] <= (m - 1) * tau[m] for m in range(2, N))
print(f"N={N}: terms below N = {len(terms)}; first {len(DATA)} match DATA: {ok_data}")
print(f"primes p != 3 below N that are not terms: {missing[:10]} (count {len(missing)})")
print(f"3 is a term: {3 in term_set}")
print(f"phi(m)/tau(m) <= (m-1)/2 for all 2 <= m < N: {ok_ineq}")
ok = ok_data and not missing and 3 not in term_set and ok_ineq
print("RESULT", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
