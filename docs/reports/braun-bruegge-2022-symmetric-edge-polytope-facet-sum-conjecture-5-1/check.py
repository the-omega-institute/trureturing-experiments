"""Braun-Bruegge, arXiv:2201.13303v4, Conjecture 5.1: exact checks of all four assertions for n <= 200.

F(a, b, c) = sum_{j=0}^{c} C(c, j) C(b, (b - c)/2 + j) C(a, (a - c)/2 + j) for a >= b >= c >= 1 of a common parity
(Proposition 4.4 with t = 3). Checks, for every admissible triple with a + b + c = n + 1 and 4 <= n <= 200:
the maximum is F(n - 1, 1, 1) for even n and F(n - 3, 2, 2) for odd n >= 5; F(a, b, c) <= F(a + 2, b, c - 2) when c >= 3;
F(a, b, c) <= F(a + 2, b - 2, c) when c >= 3 and b - 2 >= c. Also checks the sign pattern of the search seat's transfer
lemma: Delta(t) = p_{a+2}(t) p_{c-2}(t) - p_a(t) p_c(t) is >= 0 for small |t| and <= 0 for large |t|, with total sum 0.
"""
import sys
from math import comb


def F(a, b, c):
    return sum(comb(c, j) * comb(b, (b - c) // 2 + j) * comb(a, (a - c) // 2 + j) for j in range(c + 1))


def p(m, t):
    if abs(t) > m or (m + t) % 2:
        return 0
    return comb(m, (m + t) // 2)


ok = True
triples = 0
for n in range(4, 201):
    s = n + 1
    best = None
    for c in range(1, s):
        for b in range(c, s):
            a = s - b - c
            if a < b or (a - b) % 2 or (b - c) % 2:
                continue
            triples += 1
            v = F(a, b, c)
            best = v if best is None else max(best, v)
            if c >= 3 and F(a + 2, b, c - 2) < v:
                ok = False
                print("first transfer fails", a, b, c)
            if c >= 3 and b - 2 >= c and F(a + 2, b - 2, c) < v:
                ok = False
                print("second transfer fails", a, b, c)
    target = F(n - 1, 1, 1) if n % 2 == 0 else (F(n - 3, 2, 2) if n >= 5 else None)
    if target is not None and best != target:
        ok = False
        print("maximum fails at n =", n)
sign_ok = True
for a in range(3, 60):
    for c in range(3, a + 1):
        if (a - c) % 2:
            continue
        D = [p(a + 2, t) * p(c - 2, t) - p(a, t) * p(c, t) for t in range(-a - 2, a + 3)]
        if sum(D) != 0:
            sign_ok = False
        nz = [x for x in D if x != 0]
        mid = len(D) // 2
        right = [x for x in D[mid:] if x != 0]
        changes = sum(1 for u, w in zip(right, right[1:]) if (u > 0) != (w > 0))
        if changes > 1 or (right and right[-1] > 0):
            sign_ok = False
ok = ok and sign_ok
print(f"admissible triples checked: {triples} (4 <= n <= 200)")
print(f"maxima and both transfers: {ok}; transfer-lemma sign pattern for 3 <= c <= a < 60: {sign_ok}")
print(f"F(9,1,1) = {F(9, 1, 1)}, F(8,2,2) = {F(8, 2, 2)}")
print("RESULT PASS" if ok else "RESULT FAIL")
sys.exit(0 if ok else 1)
