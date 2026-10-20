"""OEIS A400584 (Mendonca): T(n,k) = number of +/min labelings of the internal vertices of the complete plane
binary tree of depth n, leaves 1, with root value k.  Conjecture: T(n,k) > T(n,k+1) for n >= 3, 1 <= k < 2^n.

1. Brute force over all labelings for n <= 4 and comparison with the recurrence
   T(n+1,k) = T(n,k)^2 + 2 T(n,k) sum_{u>k} T(n,u) + sum_{i<k} T(n,i) T(n,k-i).
2. The recurrence up to n = 12: row sums 2^(2^n - 1) and strict decrease for 3 <= n <= 12.
3. The adjacent-difference identity
   b_k - b_{k+1} = sum_{i<k} (a_i - a_k)(a_{k-i} - a_{k+1-i}) + a_{k+1}(2 a_k - a_{k+1}) + 2 (a_k - a_{k+1}) sum_{j>=k+2} a_j
   on random strictly decreasing positive sequences (fixed seed) and on rows 1..10, for every 1 <= k <= 2m.
4. The three strict cases (k < m, k = m, m < k < 2m) of the search seat's proof on the same rows.
"""
import random
from collections import Counter
from itertools import product

def brute(n):
    """Counter of root values over all labelings of depth n (internal vertices in heap order)."""
    m = 2 ** n - 1
    c = Counter()
    for lab in product((0, 1), repeat=m):
        vals = [1] * (2 ** n)
        off = m
        # evaluate bottom-up: level n-1 internal vertices are heap indices 2^(n-1)-1 .. 2^n-2
        level = vals
        for d in range(n - 1, -1, -1):
            start = 2 ** d - 1
            nxt = []
            for j in range(2 ** d):
                op = lab[start + j]
                a, b = level[2 * j], level[2 * j + 1]
                nxt.append(a + b if op else min(a, b))
            level = nxt
        c[level[0]] += 1
    return c

def step(a):
    m = len(a)
    A = [0] + list(a) + [0] * (m + 2)
    out = []
    for k in range(1, 2 * m + 1):
        s = sum(A[j] for j in range(k + 1, m + 1))
        c = sum(A[i] * A[k - i] for i in range(1, k))
        out.append(A[k] * A[k] + 2 * A[k] * s + c)
    return out

def ident_ok(a):
    m = len(a)
    A = [0] + list(a) + [0] * (2 * m + 3)
    B = [0] + step(a) + [0] * 3
    for k in range(1, 2 * m + 1):
        rhs = (sum((A[i] - A[k]) * (A[k - i] - A[k + 1 - i]) for i in range(1, k))
               + A[k + 1] * (2 * A[k] - A[k + 1])
               + 2 * (A[k] - A[k + 1]) * sum(A[j] for j in range(k + 2, m + 1)))
        if B[k] - B[k + 1] != rhs:
            return False
    return True

def cases_ok(a):
    m = len(a)
    A = [0] + list(a) + [0] * (2 * m + 3)
    for k in range(1, 2 * m):
        if k < m:
            t = A[k + 1] * (2 * A[k] - A[k + 1])
        elif k == m:
            t = (A[1] - A[m]) * (A[m - 1] - A[m])
        else:
            i = k - m
            t = (A[i] - A[k]) * (A[k - i] - A[k + 1 - i])
        if t <= 0:
            return False
    return True

ok = True
rows = [[1]]
for n in range(1, 13):
    rows.append(step(rows[-1]))
for n in range(0, 5):
    c = brute(n)
    bt = [c[k] for k in range(1, 2 ** n + 1)]
    print(f"brute n={n}: {bt[:8]}{' ...' if len(bt) > 8 else ''} matches recurrence: {bt == rows[n]}")
    ok &= bt == rows[n] and sum(c.values()) == 2 ** (2 ** n - 1) and set(c) == set(range(1, 2 ** n + 1))
for n in range(0, 13):
    r = rows[n]
    ok &= sum(r) == 2 ** (2 ** n - 1) and len(r) == 2 ** n and all(x > 0 for x in r)
    if n >= 3:
        ok &= all(r[k] > r[k + 1] for k in range(len(r) - 1))
print("rows 0..2:", rows[0], rows[1], rows[2], "row 3:", rows[3])
print("strict decrease 3 <= n <= 12:", all(all(rows[n][k] > rows[n][k + 1] for k in range(2 ** n - 1)) for n in range(3, 13)))
random.seed(400584)
rand_ok = True
for _ in range(4000):
    m = random.randint(2, 10)
    a = sorted(random.sample(range(1, 200), m), reverse=True)
    rand_ok &= ident_ok(a) and cases_ok(a)
row_ok = all(ident_ok(rows[n]) for n in range(1, 11)) and all(cases_ok(rows[n]) for n in range(3, 11))
print("identity and strict cases on 4000 random sequences:", rand_ok, "on rows:", row_ok)
ok &= rand_ok and row_ok
print("RESULT", "PASS" if ok else "FAIL")
raise SystemExit(0 if ok else 1)
