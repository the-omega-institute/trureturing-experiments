"""Binary-outcome detection loophole: the symmetric threshold conjecture eta*(m) = 4/(m+4).

Source: T. Cope, "The Binary-Outcome Detection Loophole", arXiv:2005.03344v3; New J. Phys. 23 (2021) 073032,
Section 4.3: "our conjecture of eta* = 4/(m + 4) for m_A = m_B = m".

Checks:
 1. Three-input leader model at m = 13, eta = 6/25: exact enumeration (Fractions) of all
    2 * C(13,3) * 8 = 4576 active atoms plus the all-F atom, for six Boolean matrices G
    (structured and random); the induced distribution equals D_eta p_G entry by entry
    (conclusive pairs eta^2 p_G, single failures eta(1-eta)/2, FF (1-eta)^2).
 2. The same identity at m = 14, 16, 20 for two random G each.
 3. 6/(m+12) > 4/(m+4) exactly for 13 <= m <= 200, and the gap 2(m-12)/((m+12)(m+4)).
The vertex structure of the binary no-signalling polytope (local relabelings of parity boxes) is cited from
the literature (Jones-Masanes; Barrett et al. 2005) and is not checked here.
"""
import itertools, random, sys
from fractions import Fraction as Fr

bad = 0
def ok(c, msg):
    global bad
    print(("ok   " if c else "FAIL ") + msg)
    if not c: bad += 1

def model_distribution(m, G):
    eta = Fr(6, m + 12); t = Fr(48 * m, (m + 12) ** 2); beta = 1 - t
    subsets = list(itertools.combinations(range(m), 3))
    w = t / 2 / len(subsets) / 8          # leader choice, subset, bits
    P = {}
    def add(x, y, a, b, val):
        P[(x, y, a, b)] = P.get((x, y, a, b), 0) + val
    for leader in (0, 1):
        for S in subsets:
            for bits in itertools.product((0, 1), repeat=3):
                r = dict(zip(S, bits))
                for x in range(m):
                    for y in range(m):
                        if leader == 0:
                            a = r[x] if x in S else 'F'
                            vals = {r[s] ^ G[s][y] for s in S}
                            b = vals.pop() if len(vals) == 1 else 'F'
                        else:
                            b = r[y] if y in S else 'F'
                            vals = {r[s] ^ G[x][s] for s in S}
                            a = vals.pop() if len(vals) == 1 else 'F'
                        add(x, y, a, b, w)
    for x in range(m):
        for y in range(m):
            add(x, y, 'F', 'F', beta)
    return eta, P

def target(m, G, eta):
    T = {}
    for x in range(m):
        for y in range(m):
            for a in (0, 1):
                for b in (0, 1):
                    T[(x, y, a, b)] = eta * eta * (Fr(1, 2) if a ^ b == G[x][y] else 0)
                T[(x, y, a, 'F')] = eta * (1 - eta) / 2
                T[(x, y, 'F', a)] = eta * (1 - eta) / 2
            T[(x, y, 'F', 'F')] = (1 - eta) ** 2
    return T

def same(m, G):
    eta, P = model_distribution(m, G)
    T = target(m, G, eta)
    keys = set(P) | set(T)
    return all(P.get(k, 0) == T.get(k, 0) for k in keys)

rng = random.Random(20261012)
m = 13
Gs = {
    "popcount(x&y) mod 2": [[bin(x & y).count('1') % 2 for y in range(m)] for x in range(m)],
    "x==y": [[1 if x == y else 0 for y in range(m)] for x in range(m)],
    "quadratic mod 13 < 6": [[1 if (x * x + 3 * x * y + 5 * y * y + 2 * x + y) % 13 < 6 else 0 for y in range(m)] for x in range(m)],
    "all zero": [[0] * m for _ in range(m)],
}
for k in range(2):
    Gs[f"random {k}"] = [[rng.randint(0, 1) for _ in range(m)] for _ in range(m)]
for name, G in Gs.items():
    ok(same(m, G), f"m=13, G={name}: local model equals D_eta p_G exactly (eta=6/25)")
for m2 in (14, 16, 20):
    for k in range(2):
        G = [[rng.randint(0, 1) for _ in range(m2)] for _ in range(m2)]
        ok(same(m2, G), f"m={m2}, random G {k}: local model equals D_eta p_G exactly")
ok(all(Fr(6, mm + 12) > Fr(4, mm + 4) and Fr(6, mm + 12) - Fr(4, mm + 4) == Fr(2 * (mm - 12), (mm + 12) * (mm + 4)) for mm in range(13, 201)),
   "6/(m+12) > 4/(m+4) with gap 2(m-12)/((m+12)(m+4)) for 13 <= m <= 200")
ok(Fr(6, 12 + 12) == Fr(4, 12 + 4) and Fr(6, 12 + 12) == Fr(1, 4), "equality at m = 12 (both 1/4)")

print("bad=", bad)
sys.exit(1 if bad else 0)
