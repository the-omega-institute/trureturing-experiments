"""OEIS A375786 conjecture: among integer cuboids x <= y <= z with surface 6n^2
(xy + yz + zx = 3n^2), the one that is lexicographically least in (x, y) has the
minimum volume.

For every n <= N the script enumerates all solutions via divisors of x^2 + 3n^2,
checks the conjecture, compares the minimum with the published data, and checks
the lemmas of the proof: 2xyz >= x(S - x^2) for every solution; the least edge
x0 is 1 or 2, and x0 = 2 only when 1 + 3n^2 is prime, n is even and y0 = 2.
"""
import sys
from math import isqrt

N = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
DATA = [1,8,13,36,37,104,73,188,121,252,181,428,253,540,337,764,433,828,541,1196,661,1448,
        793,1476,937,2024,1093,2160,1261,2592,1441,2628,1633,3464,1837,3884]

def is_prime(m):
    return m >= 2 and all(m % q for q in range(2, isqrt(m) + 1))

def solutions(n):
    S = 3 * n * n
    out = []
    for x in range(1, n + 1):          # 3x^2 <= S
        T = x * x + S
        for d in range(2 * x, isqrt(T) + 1):
            if T % d == 0:
                out.append((x, d - x, T // d - x))
    return out

ok_conj = ok_data = ok_lemma = ok_x0 = True
for n in range(1, N + 1):
    S = 3 * n * n
    sol = solutions(n)
    assert all(x <= y <= z and x * y + y * z + z * x == S for x, y, z in sol)
    lex = min(sol, key=lambda t: (t[0], t[1]))
    vmin = min(x * y * z for x, y, z in sol)
    if lex[0] * lex[1] * lex[2] != vmin:
        ok_conj = False
        print("counterexample", n, lex, vmin)
    if n <= len(DATA) and vmin != DATA[n - 1]:
        ok_data = False
    if any(2 * x * y * z < x * (S - x * x) for x, y, z in sol):
        ok_lemma = False
    x0 = lex[0]
    if x0 == 2:
        if not (is_prime(1 + S) and n % 2 == 0 and lex[1] == 2):
            ok_x0 = False
    elif x0 != 1:
        ok_x0 = False
print(f"N={N}: conjecture holds for all n <= N: {ok_conj}")
print(f"minimum volume matches the first {len(DATA)} terms: {ok_data}")
print(f"2xyz >= x(S - x^2) on every solution: {ok_lemma}")
print(f"least edge is 1, or 2 with 1 + 3n^2 prime, n even, y0 = 2: {ok_x0}")
good = ok_conj and ok_data and ok_lemma and ok_x0
print("RESULT", "PASS" if good else "FAIL")
sys.exit(0 if good else 1)
