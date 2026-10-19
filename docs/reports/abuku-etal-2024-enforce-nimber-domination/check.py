"""Checks for Abuku-Kimura-Kiya-Larsson-Saha-Suetsugu-Yamashita, arXiv:2311.01006v2, Problem 42: domination A |- B
(O_{A.B} = O_A everywhere) does not imply equality of the enforce nimbers with the nimbers of A.
Outcomes and nimbers are computed by the recursions of Definitions 5, 13, 29 and 33 on finite acyclic boards."""
from itertools import product, combinations

def mex(s):
    k = 0
    while k in s:
        k += 1
    return k

def solve(n, fA, fB):
    GA, GE, OA, OE = {}, {}, {}, {}
    for x in range(n):  # options always go to smaller positions, so this order is a valid recursion order
        GA[x] = mex({GA[y] for y in fA[x]})
        GE[x] = mex({GE[y] for y in fA[x]} & {GE[y] for y in fB[x]})
        OA[x] = 'P' if all(OA[y] == 'N' for y in fA[x]) else 'N'
        OE[x] = 'P' if (all(OE[y] == 'N' for y in fA[x]) or all(OE[y] == 'N' for y in fB[x])) else 'N'
    return GA, GE, OA, OE

ok = True
fA = [[], [0], [0, 1]]; fB = [[], [0], [0]]
GA, GE, OA, OE = solve(3, fA, fB)
w = OA == OE and GA[2] == 2 and GE[2] == 1
print("three-position witness: OA", [OA[x] for x in range(3)], "OE", [OE[x] for x in range(3)],
      "GA", [GA[x] for x in range(3)], "GE", [GE[x] for x in range(3)], w)
ok &= w

def family(d, N):
    fA = [[m for m in range(N + 1) if m + 2 <= n] for n in range(N + 1)]
    fB = [[d * k for k in range(N + 1) if d * k < n] for n in range(N + 1)]
    return solve(N + 1, fA, fB)
fam = True
for d in range(3, 8):
    for N in range(2 * d, 60):
        GA, GE, OA, OE = family(d, N)
        fam &= OA == OE and all(GA[n] == n // 2 for n in range(N + 1))
        fam &= all(GE[n] == (0 if n < 2 else -(-n // d)) for n in range(N + 1))
        fam &= all(GA[2 * d * r] - GE[2 * d * r] == (d - 2) * r for r in range(1, N // (2 * d) + 1))
print("family (subtract >= 2 vs. multiples of d), d = 3..7, N = 2d..59: domination, closed forms, gap (d-2)r:", fam)
ok &= fam

# minimality: every short ruleset on at most two positions is acyclic, hence (after relabelling) has options to smaller positions only
def maps(n):
    return [[list(c) for r in range(x + 1) for c in combinations(range(x), r)] for x in range(n)]
small = True
for n in (1, 2):
    for a in product(*maps(n)):
        for b in product(*maps(n)):
            GA, GE, OA, OE = solve(n, list(a), list(b))
            if OA == OE and GA != GE:
                small = False
print("no counterexample on boards with at most two positions:", small)
ok &= small
print("RESULT", "PASS" if ok else "FAIL")
