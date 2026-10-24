"""Sign of the Turanian Delta_f(x) on x = 0, 1/4, ..., 99.75 for p = 0, b = (B,...,B), mu = 1, alpha = beta = 1/2.
Usage: python3 turanian_scan.py B Q [Q ...]"""
import sys
import mpmath as mp
from common import f, turanian

mp.mp.dps = 50
B = mp.mpf(sys.argv[1])
for q in map(int, sys.argv[2:]):
    b = [B] * q
    best = min((turanian([], b, 1, mp.mpf(1) / 2, mp.mpf(1) / 2, mp.mpf(i) / 4, 400) / f([], b, 1, mp.mpf(i) / 4, 400) ** 2, i / 4) for i in range(400))
    print(f"B={sys.argv[1]} q={q} min Delta/f^2={mp.nstr(best[0], 6)} at x={best[1]}")
