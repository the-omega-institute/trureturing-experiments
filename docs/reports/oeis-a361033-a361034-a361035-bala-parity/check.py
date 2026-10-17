"""OEIS A361033 / A361034 / A361035 (P. Bala): a_j(n) = C_j * (4n)!/(n! (n+j)!^3) with
C_1 = 3, C_2 = 2520, C_3 = 9979200. Conjectures: a_j(n) is odd iff n + j is a power of 2
(with n = 2^k - j, k >= 0, 1, 2 respectively).

For n <= N the script checks that the quotients are integers, the parity conjectures, the
published data, and the 2-adic valuation v2(a_j(n)) = 3*(s2(n + j) - 1) used in the proof
(s2 = binary digit sum).
"""
import sys
from math import factorial

N = int(sys.argv[1]) if len(sys.argv) > 1 else 600
C = {1: 3, 2: 2520, 3: 9979200}
DATA = {1: [3, 9, 280, 17325, 1513512, 162954792, 20193091776, 2768662192725],
        2: [315, 280, 3675, 116424, 5885880, 399072960, 33129291195, 3190228041000],
        3: [46200, 17325, 116424, 2134440, 67953600, 3086579925, 179961581800]}

def v2(x):
    return (x & -x).bit_length() - 1

ok = True
for j in (1, 2, 3):
    vals = []
    for n in range(N + 1):
        num = C[j] * factorial(4 * n)
        den = factorial(n) * factorial(n + j) ** 3
        if num % den: ok = False; print("not integral", j, n); break
        a = num // den
        vals.append(a)
        pow2 = (n + j) & (n + j - 1) == 0
        if (a % 2 == 1) != pow2: ok = False; print("parity", j, n)
        if v2(a) != 3 * (bin(n + j).count("1") - 1): ok = False; print("v2", j, n)
    if vals[:len(DATA[j])] != DATA[j]: ok = False; print("data", j)
print(f"N={N}: integral, parity iff n+j is a power of 2, v2 = 3(s2(n+j)-1), data match, for j = 1, 2, 3: {ok}")
print("RESULT", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
