"""Random control: p <= q <= 4 with condition (13), mu >= 1, alpha, beta in [0,3], m <= 15.
Usage: python3 random_control.py SEED DRAWS"""
import random
import sys
import mpmath as mp
from common import delta, cond13, c

mp.mp.dps = 40
seed, draws = int(sys.argv[1]), int(sys.argv[2])
random.seed(seed)
accepted = violations = 0
for _ in range(draws):
    p = random.randint(0, 3)
    q = random.randint(max(p, 1), 4)
    a = [mp.mpf(random.choice([random.uniform(0.01, 1), random.uniform(0.01, 6)])) for _ in range(p)]
    b = [mp.mpf(random.choice([random.uniform(0.01, 1), random.uniform(0.01, 6)])) for _ in range(q)]
    if not cond13(a, b):
        continue
    accepted += 1
    mu = mp.mpf(random.choice([1, random.uniform(1, 1.3), random.uniform(1, 4)]))
    al, be = mp.mpf(random.uniform(0, 3)), mp.mpf(random.uniform(0, 3))
    scale = abs(c(a, b, mu) * c(a, b, mu + al + be))
    for m in range(16):
        if delta(a, b, mu, al, be, m) / (scale * 2 ** m / mp.factorial(m)) < -mp.mpf("1e-25"):
            violations += 1
            print("violation", [float(x) for x in a], [float(x) for x in b], float(mu), float(al), float(be), m)
            break
print(f"seed={seed} draws={draws} accepted={accepted} violations={violations}")
