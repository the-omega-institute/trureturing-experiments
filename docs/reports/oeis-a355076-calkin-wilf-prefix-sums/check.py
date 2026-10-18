"""OEIS A355076: a(n) is the denominator of S(n) = Sum_{k=0..n} fusc(k)/fusc(k+1) (fusc = Stern's diatomic A002487).
Conjecture (R. Sigrist): a(n) = 1 for n of the form 2*4^k - 1 or 2*4^k - 2.

Proof checked here: fusc(k)/fusc(k+1) runs through the Calkin-Wilf tree by rows; the row sums are
R_m = sum_{k=2^m}^{2^(m+1)-1} fusc(k)/fusc(k+1) = 3*2^(m-1) - 1/2 (children a/(a+b), (a+b)/b and the reciprocal
symmetry of each row), so S(2^(m+1) - 1) = (3*2^(m+1) - m - 4)/2, an integer iff m is even; and the last term of row m
is fusc(2^(m+1)-1)/fusc(2^(m+1)) = m + 1.
The script checks the row sums, the reciprocal symmetry, the prefix formula and the conjecture for 2^(m+1) <= N.
"""
import sys
from fractions import Fraction
N = int(sys.argv[1]) if len(sys.argv) > 1 else 1 << 16
f = [0, 1]
for k in range(2, N + 3):
    f.append(f[k // 2] if k % 2 == 0 else f[k // 2] + f[k // 2 + 1])
DATA = [1, 1, 2, 2, 6, 3, 1, 1, 4, 12, 60, 60, 12, 4, 2, 2, 10, 20, 140, 420, 840, 840, 840, 840]
S = Fraction(0); den = []
pref = {}
for n in range(0, N + 1):
    S += Fraction(f[n], f[n + 1]); den.append(S.denominator); pref[n] = S
ok_data = den[:len(DATA)] == DATA
ok_rows = True; ok_sym = True; ok_pref = True; ok_last = True
m = 0
while (1 << (m + 1)) - 1 <= N:
    row = [Fraction(f[k], f[k + 1]) for k in range(1 << m, 1 << (m + 1))]
    if sum(row) != Fraction(3 * (1 << m), 2) - Fraction(1, 2): ok_rows = False
    if sorted(row) != sorted(1 / q for q in row): ok_sym = False
    if pref[(1 << (m + 1)) - 1] != Fraction(3 * (1 << (m + 1)) - m - 4, 2): ok_pref = False
    if Fraction(f[(1 << (m + 1)) - 1], f[1 << (m + 1)]) != m + 1: ok_last = False
    m += 1
conj = all(den[2 * 4 ** k - 1] == 1 and den[2 * 4 ** k - 2] == 1 for k in range(0, 9) if 2 * 4 ** k - 1 <= N)
print(f"data: {ok_data}; row sums 3*2^(m-1)-1/2: {ok_rows}; reciprocal symmetry: {ok_sym}; prefix formula: {ok_pref}; last term m+1: {ok_last}; conjecture for 2*4^k-1 <= {N}: {conj}")
good = ok_data and ok_rows and ok_sym and ok_pref and ok_last and conj
print("RESULT", "PASS" if good else "FAIL")
sys.exit(0 if good else 1)
