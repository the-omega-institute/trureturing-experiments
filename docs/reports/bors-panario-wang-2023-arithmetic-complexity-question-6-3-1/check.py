"""Bors-Panario-Wang (arXiv:2304.00181v1) Question 6.3.1: lim inf AC(T_m)/(log^2 m/log log m) > 0.
AC(T_m) is the least number of m-congruences x = b (mod a), a | m, whose truth vectors separate all m residues.
Checks (exact search):
  1. Lemma 1 of the proposed proof, AC(T_m) >= sum over prime powers q || m of phi(q): no separating family of size
     sum phi - 1 exists, for every m <= 30 where the search is within 3e6 families;
  2. exact AC(T_m) where the search finishes, compared with the source's upper bound sum (q - 1);
  3. along primorials m = p_K#, the ratio sum phi(q) / (log^2 m / log log m) for K <= 200 (stays near 1/2).
"""
from itertools import combinations
from math import comb, log
from sympy import factorint, totient, prime

def congruences(m):
    out = set()
    for a in range(1, m + 1):
        if m % a == 0:
            for b in range(a):
                v = tuple(1 if x % a == b else 0 for x in range(m))
                if 0 < sum(v) < m:
                    out.add(v)
    return sorted(out)

def separates(fam, m):
    return len({tuple(v[x] for v in fam) for x in range(m)}) == m

def exists_family(m, K, C):
    return any(separates(f, m) for f in combinations(C, K))

ok = True
for m in range(2, 31):
    f = factorint(m)
    low = sum(int(totient(p ** a)) for p, a in f.items())
    up = sum(p ** a - 1 for p, a in f.items())
    C = congruences(m)
    if low - 1 >= 1 and comb(len(C), low - 1) > 3_000_000:
        print(f"m={m}: sum_phi={low} upper={up} (lemma search skipped: {comb(len(C), low - 1)} families)")
        continue
    lemma = (low - 1 < 1) or not exists_family(m, low - 1, C)
    ok &= lemma
    ac = None
    for K in range(low, up + 1):
        if comb(len(C), K) > 3_000_000:
            break
        if exists_family(m, K, C):
            ac = K
            break
    print(f"m={m}: sum_phi={low} upper={up} no_family_of_size_sum_phi-1={lemma} AC={ac}"
          + (f" AC==upper={ac == up}" if ac is not None else ""))
ratios = []
m_log = 0.0; s = 0
for K in range(1, 201):
    p = prime(K); m_log += log(p); s += p - 1
    if K >= 5:
        ratios.append((K, s / (m_log ** 2 / log(m_log))))
print("primorial ratios (K, sum phi / (log^2 m / log log m)) at K = 5, 50, 100, 200:",
      [(k, round(r, 4)) for k, r in ratios if k in (5, 50, 100, 200)])
ok &= min(r for _, r in ratios) > 0.3
print("RESULT", "PASS" if ok else "FAIL")
raise SystemExit(0 if ok else 1)
