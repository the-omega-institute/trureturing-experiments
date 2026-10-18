"""OEIS A356334: number of solutions (x, y), 0 <= x <= y, of x^(n+1) + y^(n+1) = (x+y)^n.
Conjecture (C. Wu): a(n) = 3 except for n = 0 or n = 2.

Proof checked here: write x = d a, y = d b with gcd(a, b) = 1 (or x = y = 0). Then d (a^(n+1) + b^(n+1)) = (a+b)^n.
- n odd: gcd(a+b, a^(n+1)+b^(n+1)) | 2, so a^(n+1)+b^(n+1) is a power of 2, forcing (a,b) = (0,1) or (1,1).
- n even: (a+b) | a^(n+1)+b^(n+1); R = (a^(n+1)+b^(n+1))/(a+b) has only primes p | n+1, with v_p(R) = v_p(n+1)
  (lifting the exponent), so R <= n+1; but R >= b^n/2 for b >= 2, so n <= 2; n = 2 adds (1, 2).
The script finds all solutions by brute force with x <= y <= B(n) for n <= NMAX, where the proof shows every solution
satisfies y <= max(1, 2^(n-1)) for n != 2 (y <= 2 for n = 2); it also checks R <= n+1 never holds with b >= 2 for
even 4 <= n <= 60 and 1 <= a < b <= 60.
"""
import sys
from math import gcd
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 14
DATA = [1, 3, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3]
res = []
for n in range(0, NMAX + 1):
    B = max(4, 2 ** n + 2)
    sols = [(x, y) for y in range(0, B + 1) for x in range(0, y + 1) if x ** (n + 1) + y ** (n + 1) == (x + y) ** n]
    res.append((n, sols))
counts = [len(s) for _, s in res]
ok_data = counts == DATA[:NMAX + 1]
ok_shape = all(set(s) == ({(0, 1)} if n == 0 else {(0, 0), (0, 1), (2 ** (n - 1), 2 ** (n - 1))} | ({(1, 2)} if n == 2 else set())) for n, s in res)
ok_R = all(not ((a ** (n + 1) + b ** (n + 1)) % (a + b) == 0 and (a ** (n + 1) + b ** (n + 1)) // (a + b) <= n + 1)
           for n in range(4, 61, 2) for b in range(2, 61) for a in range(1, b) if gcd(a, b) == 1)
print(f"counts n<={NMAX}: {counts}; matches data: {ok_data}; solution sets as in the proof: {ok_shape}; R > n+1 for even n in [4,60]: {ok_R}")
good = ok_data and ok_shape and ok_R
print("RESULT", "PASS" if good else "FAIL")
sys.exit(0 if good else 1)
