"""Ellis-Thanatipanonda, arXiv:2203.02457v1 (DAM 322 (2022)), Section 7, Family B, Conjecture 2:
for C = {1} u X u Y (X nonempty even >= 4, Y nonempty odd >= 5, x = min X, y = min Y), if 3x < y then
G_C(n) = G_{1,x}(n) for n >= 1.

The script computes nim-sequences of cut by a bitset composition DP. It checks the source's Family D list
(C = {1, 2}, 36 terms), then compares G_C with G_{1,x} for n <= N on 22 Family B cut-sets with 3x < y (fixed
three-element sets and random larger X, Y), and checks that the threshold is sharp: y = 3x - 1 first differs at
n = 4x (x = 4: n = 16; x = 6: n = 24).
"""
import numpy as np, itertools, random, sys
M = 512
idx = np.arange(M)

def nimseq(C, N):
    C = sorted(C); P = max(C) + 1
    G = [0] * (N + 1)
    # R[p][m]: bool array of xor values of compositions of m into p positive parts
    R = [[None] * (N + 1) for _ in range(P + 1)]
    for n in range(1, N + 1):
        # compute R[p][n] for p >= 2 using G[1..n-1]
        for p in range(2, P + 1):
            a = np.zeros(M, dtype=bool)
            for h in range(1, n - p + 2):
                prev = R[p - 1][n - h]
                if prev is not None:
                    a |= prev[idx ^ G[h]]
            R[p][n] = a
        opts = np.zeros(M, dtype=bool)
        for d in C:
            if d + 1 <= n and R[d + 1][n] is not None:
                opts |= R[d + 1][n]
        g = int(np.argmin(opts)) if not opts.all() else M
        G[n] = g
        a = np.zeros(M, dtype=bool); a[g] = True
        R[1][n] = a
    return G

N = int(sys.argv[1]) if len(sys.argv) > 1 else 80
ok = True
src = [0, 1, 2, 3, 1, 4, 3, 2, 4, 5, 6, 7, 8, 9, 7, 6, 9, 8, 11, 10, 12, 13, 10, 11, 13, 12, 15, 14, 16, 17, 5, 15, 17,
       16, 19, 18]
fd = nimseq({1, 2}, 36)
print("Family D (C = {1, 2}) matches the source's 36 terms:", fd[1:37] == src)
ok &= fd[1:37] == src
random.seed(1)
cases = []
for x in (4, 6, 8):
    for y0 in range(3 * x + 1, 3 * x + 8):
        if y0 % 2 == 1:
            cases.append((x, {1, x, y0}))
for _ in range(10):
    x = random.choice((4, 6))
    Xe = {x} | {e for e in range(x + 2, 40, 2) if random.random() < 0.2}
    y = random.choice([v for v in range(3 * x + 1, 3 * x + 12) if v % 2 == 1])
    Yo = {y} | {o for o in range(y + 2, 45, 2) if random.random() < 0.2}
    cases.append((x, {1} | Xe | Yo))
base = {}
bad = 0
for x, C in cases:
    if x not in base:
        base[x] = nimseq({1, x}, N)
    G = nimseq(C, N)
    diff = [n for n in range(1, N + 1) if G[n] != base[x][n]]
    print(f"x={x} C={sorted(C)}: first diff {diff[:3]}")
    bad += bool(diff)
print("Family B cut-sets with 3x < y:", len(cases), "violations up to n =", N, ":", bad)
ok &= bad == 0 and len(cases) == 22
for x in (4, 6):
    y = 3 * x - 1
    G = nimseq({1, x, y}, 4 * x + 4)
    b = nimseq({1, x}, 4 * x + 4)
    d = [n for n in range(1, 4 * x + 5) if G[n] != b[n]]
    print(f"boundary x={x} y={y}: first diff n={d[0] if d else None} (expected {4 * x})")
    ok &= bool(d) and d[0] == 4 * x
print("RESULT PASS" if ok else "RESULT FAIL")
sys.exit(0 if ok else 1)
