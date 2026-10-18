"""Bala's fractional factorial ratios, OEIS A347854-A347858 ("conjectured to be integral"), and the b(n), c(n)
of A347855 ("appears to be integral"). x! := Gamma(1 + x).

Exact arithmetic only:
  Gamma(1 + t/2)/sqrt(pi) = t!/(2^t ((t-1)/2)!) for odd t (sqrt(pi) cancels in every ratio);
  Gamma(1 + 4n/3)/Gamma(1 + n/3) = prod_{j=1..n} (n/3 + j).
Checks: the published data, integrality for n <= N, the relations a(3m+1) = 3^(3m+1) b(m) and
a(3m+2) = 3^(3m+2) c(m), and nonnegativity of the ten Landau-type step functions on a common grid of period 1.
Usage: python3 check.py N
"""
import sys
from fractions import Fraction
from math import factorial as F, floor, prod

def g_half(t):
    return Fraction(F(t // 2)) if t % 2 == 0 else Fraction(F(t), 2 ** t * F((t - 1) // 2))

def a854(n): return Fraction(F(6 * n), F(3 * n) * F(2 * n)) * g_half(n) / g_half(3 * n)
def a856(n): return Fraction(F(6 * n), F(4 * n) * F(n)) * g_half(n) / g_half(3 * n)
def a857(n): return Fraction(F(6 * n), F(3 * n) * F(2 * n)) * g_half(3 * n) / g_half(5 * n)
def a858(n): return Fraction(F(9 * n), F(3 * n) * F(2 * n)) * g_half(n) / g_half(9 * n)
def a855(n):
    p = Fraction(1)
    for j in range(1, n + 1): p *= Fraction(n, 3) + j
    return Fraction(F(4 * n), F(2 * n) * F(n)) / p
def b(m): return Fraction(F(12 * m + 4), F(6 * m + 2) * F(3 * m + 1) * prod(3 * m + 3 * k + 4 for k in range(3 * m + 1)))
def c(m): return Fraction(F(12 * m + 8), F(6 * m + 4) * F(3 * m + 1) * prod(3 * m + 3 * k + 2 for k in range(3 * m + 3)))

DATA = {
    'A347854': (a854, [1, 40, 4620, 622336, 89237148, 13236695040, 2005604901300, 308350245273600]),
    'A347855': (a855, [1, 9, 189, 4620, 120285, 3241134, 89237148, 2493521172, 70429218525]),
    'A347856': (a856, [1, 20, 990, 56576, 3432198, 215147520, 13768454700, 893768826880]),
    'A347857': (a857, [1, 24, 1386, 91392, 6374082, 458625024, 33660501840, 2504739913728]),
    'A347858': (a858, [1, 512, 1021020, 2399141888, 6016814703900, 15626259253952512]),
}

def fl(x): return floor(x)
STEP = {
    'H854': lambda x: fl(12*x) + fl(x) - fl(6*x) - fl(4*x) - fl(3*x),
    'H856': lambda x: fl(12*x) + fl(x) - fl(8*x) - fl(2*x) - fl(3*x),
    'H857': lambda x: fl(12*x) + fl(3*x) - fl(6*x) - fl(4*x) - fl(5*x),
    'H858': lambda x: fl(18*x) + fl(x) - fl(6*x) - fl(4*x) - fl(9*x),
    'K854': lambda x: fl(12*x) + fl(2*x) + fl(3*x) - 2*fl(6*x) - fl(4*x) - fl(x),
    'K856': lambda x: fl(12*x) + fl(3*x) - fl(8*x) - fl(6*x) - fl(x),
    'K857': lambda x: fl(12*x) + fl(5*x) - fl(4*x) - fl(10*x) - fl(3*x),
    'K858': lambda x: fl(2*x) + fl(9*x) - fl(6*x) - fl(4*x) - fl(x),
    'T1': lambda x: fl(12*x) - fl(6*x) - fl(3*x) - fl(4*x - Fraction(1, 3)) + fl(x - Fraction(1, 3)),
    'T2': lambda x: fl(12*x) - fl(6*x) - fl(3*x) - fl(4*x - Fraction(2, 3)) + fl(x - Fraction(2, 3)),
}

def main():
    N = int(sys.argv[1]); ok = True
    for name, (f, data) in DATA.items():
        d = [f(n) for n in range(len(data))] == data
        bad = [n for n in range(N + 1) if f(n).denominator != 1]
        print(f"{name}: data {d}; non-integral n <= {N}: {bad[:5]}"); ok &= d and bad == []
    rb = all(a855(3 * m + 1) == 3 ** (3 * m + 1) * b(m) for m in range(N // 3 + 1))
    rc = all(a855(3 * m + 2) == 3 ** (3 * m + 2) * c(m) for m in range(N // 3 + 1))
    ib = [m for m in range(N // 3 + 1) if b(m).denominator != 1]; ic = [m for m in range(N // 3 + 1) if c(m).denominator != 1]
    print(f"A347855 b, c: relations to a {rb and rc}; non-integral b, c for m <= {N // 3}: {ib[:5]}, {ic[:5]}"); ok &= rb and rc and ib == [] and ic == []
    L = 360  # multiple of every floor coefficient and of the shifted breakpoints (3j + r)/12
    for name, s in STEP.items():
        vals = [s(Fraction(j, L)) for j in range(L)]
        per = all(s(Fraction(j, L) + 1) == v for j, v in enumerate(vals))
        print(f"{name}: min {min(vals)} max {max(vals)} on the grid j/{L}; period 1: {per}"); ok &= min(vals) >= 0 and per
    print("RESULT", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)

main()
