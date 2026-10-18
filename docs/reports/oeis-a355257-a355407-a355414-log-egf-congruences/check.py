"""OEIS A355257 (array A(n,k) = k! [x^k] log((1-x)/(1-2x)) / (1-x)^n), A355407 (n = 4), A355414 (n = 5).

Checked statements:
- C1 (A355257, A355407, A355414): for p prime, A(n, p) = -1 (mod p).
- C2 (A355257): for n >= 0, k >= 1, k != 4: k | A(n, k) iff k is not prime.
- C4 (A355257): for p prime and n even, A(n, p)/(p - 1) is an integer = 1 (mod p).

Proof ingredients checked here:
- A(n, k) = sum_{j=1..k} (k!/j) (2^j - 1) C(k-j+n-1, n-1), from log((1-x)/(1-2x)) = sum_j (2^j-1) x^j / j;
  this is compared with a direct series expansion (sympy) for n <= 6, k <= 12.
- The j < k terms are divisible by k (j | (k-1)!), and the j = k term is (k-1)!(2^k - 1).
- Modulo p - 1 only the j = p - 1 term survives: p (p-2)! (2^(p-1) - 1) n.
"""
import sys
from math import comb, factorial
from fractions import Fraction
from sympy import isprime, symbols, log, series, Rational, primerange
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 30
KMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 200
def A(n, k):
    s = 0
    for j in range(1, k + 1):
        c = comb(k - j + n - 1, n - 1) if n > 0 else (1 if j == k else 0)
        s += (factorial(k) // j) * (2 ** j - 1) * c
    return s
x = symbols('x')
ok_series = True
for n in range(0, 7):
    ser = series(log((1 - x) / (1 - 2 * x)) / (1 - x) ** n, x, 0, 13).removeO()
    for k in range(0, 13):
        if ser.coeff(x, k) * factorial(k) != A(n, k): ok_series = False
D355407 = [0, 1, 11, 110, 1154, 13144, 164136, 2251920, 33923760, 560180160, 10117886400]
D355414 = [0, 1, 13, 149, 1750, 21894, 295500, 4320420, 68487120, 1176564240, 21883528800]
D355257 = [0, 0, 1, 0, 1, 3, 0, 1, 5, 14, 0, 1, 7, 29, 90, 0, 1, 9, 50, 206, 744, 0, 1, 11, 77, 406, 1774, 7560]
anti = [A(d - k, k) for d in range(0, 7) for k in range(0, d + 1)]
ok_data = [A(4, k) for k in range(11)] == D355407 and [A(5, k) for k in range(11)] == D355414 and anti == D355257
P = list(primerange(2, KMAX + 1))
c1 = all(A(n, p) % p == p - 1 for n in range(0, NMAX + 1) for p in P)
c2 = all(((A(n, k) % k == 0) == (not isprime(k))) for n in range(0, 16) for k in range(1, KMAX + 1) if k != 4)
c4 = all(A(n, p) % (p - 1) == 0 and (A(n, p) // (p - 1)) % p == 1 for n in range(0, NMAX + 1, 2) for p in P)
mod_pm1 = all((A(n, p) - p * factorial(p - 2) * (2 ** (p - 1) - 1) * n) % (p - 1) == 0 for n in range(0, NMAX + 1) for p in P if p > 2)
k4 = [A(n, 4) % 4 for n in range(0, 8)]
print(f"series check n<=6, k<=12: {ok_series}; data A355257/A355407/A355414: {ok_data}")
print(f"C1 n<={NMAX}, p<={KMAX}: {c1}; C2 n<16, k<={KMAX}: {c2}; C4 even n<={NMAX}: {c4}; reduction mod p-1: {mod_pm1}")
print(f"k = 4 (excluded in C2): A(n,4) mod 4 for n = 0..7: {k4}")
good = ok_series and ok_data and c1 and c2 and c4 and mod_pm1
print("RESULT", "PASS" if good else "FAIL")
sys.exit(0 if good else 1)
