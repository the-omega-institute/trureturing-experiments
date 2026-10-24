"""Conjecture 2 witness: a1 = b1 = (), a2 = (1), b2 = (1/2)^q, mu = 1, alpha = beta = 1/2.
Compares the closed form (45 pi/16) 2^-q - 4 (16/(9 pi))^q - 6 (16/(15 pi))^q with the direct
Cauchy sum of the first Turanian coefficient for q = 1..11."""
import mpmath as mp
from common import delta

mp.mp.dps = 40
h = mp.mpf(1) / 2
pi = mp.pi
worst = mp.mpf(0)
for q in range(1, 12):
    closed = 45 * pi / 16 * mp.mpf(2) ** (-q) - 4 * (16 / (9 * pi)) ** q - 6 * (16 / (15 * pi)) ** q
    direct = delta([mp.mpf(1)], [h] * q, 1, h, h, 1)
    worst = max(worst, abs(closed - direct))
    print(f"q={q} closed={mp.nstr(closed, 10)} direct={mp.nstr(direct, 10)}")
print(f"max |closed-direct| = {mp.nstr(worst, 3)}")
