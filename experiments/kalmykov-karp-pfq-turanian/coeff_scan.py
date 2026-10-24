"""First Turanian coefficient for p = 0 and repeated b_j, closed form and direct sums."""
import sys
import mpmath as mp
from common import delta

mp.mp.dps = 40
h = mp.mpf(1) / 2
print("closed form, b=1/2, mu=1, alpha=beta=1/2: delta_1 = 2*2^-q - (16/(9 pi))^q - (16/(15 pi))^q")
for q in range(1, 9):
    v = 2 * mp.mpf(2) ** (-q) - (16 / (9 * mp.pi)) ** q - (16 / (15 * mp.pi)) ** q
    d = delta([], [h] * q, 1, h, h, 1)
    print(f"q={q} closed={mp.nstr(v, 10)} direct={mp.nstr(d, 10)}")
print("least q with delta_1<0, b=1/2, alpha=beta=1/2")
for mu in [1, 1.5, 2, 3]:
    q = next(q for q in range(1, 40) if delta([], [h] * q, mu, h, h, 1) < 0)
    print(f"mu={mu} q={q}")
print("sign of delta_1, q=6, b=1/2, mu=1; rows alpha, columns beta in 0.1 0.25 0.5 0.75 1 1.5 2")
for al in [0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0]:
    row = "".join("-" if delta([], [h] * 6, 1, al, be, 1) < 0 else "+" for be in [0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0])
    print(f"alpha={al} {row}")
print("p=1, b=(1/2)^6, mu=1, alpha=beta=1/2")
for a in [1, 2, 3]:
    r = delta([mp.mpf(a)], [h] * 6, 1, h, h, 1) / delta([mp.mpf(a)], [h] * 6, 1, h, h, 0)
    print(f"a={a} delta_1/delta_0={mp.nstr(r, 6)}")
sys.exit(0)
