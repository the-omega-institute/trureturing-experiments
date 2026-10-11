#!/usr/bin/env python3
"""Checks for Gabric-Shallit, arXiv:2302.13147v2, Conjecture 5 (expected width of the
largest block-palindrome factorization).  Usage: python3 check.py B274190 B274199
(the OEIS b-files of A274190 and A274199).  Exits 0 and prints `bad= 0` when every check passes."""
import itertools, sys
from fractions import Fraction
from decimal import Decimal, getcontext

bad = 0
def ok(cond, msg):
    global bad
    if not cond:
        bad += 1
        print("FAIL", msg)

# --- definition: largest BP-factorization width, by exhaustive search -------------
def maxw_def(w):
    best = 1 if len(w) > 0 else 0
    for l in range(1, len(w) // 2 + 1):
        if w[:l] == w[len(w) - l:]:
            best = max(best, 2 + maxw_def(w[l:len(w) - l]))
    return best

def maxw_greedy(w):
    r = 0
    while True:
        for l in range(1, len(w) // 2 + 1):
            if w[:l] == w[len(w) - l:]:
                w = w[l:len(w) - l]; r += 1
                break
        else:
            return 2 * r + (1 if len(w) > 0 else 0)

for k, N in ((2, 14), (3, 9), (4, 7)):
    for n in range(0, N + 1):
        for w in itertools.product(range(k), repeat=n):
            ok(maxw_def(w) == maxw_greedy(w), f"greedy k={k} w={w}")

# --- source Table 2 (IB_2(n,t), n=10..20, t=1..10) via the greedy width -------------
table = {10: [284,12,224,40,168,72,96,64,32,32], 11: [568,0,472,0,416,0,336,0,192,0],
         12: [1116,20,856,88,656,176,448,224,224,160], 13: [2232,0,1752,0,1488,0,1248,0,896,0],
         14: [4424,40,3328,176,2544,432,1856,640,1152,640]}
for n, row in table.items():
    cnt = [0] * (n + 1)
    for w in itertools.product(range(2), repeat=n):
        cnt[maxw_greedy(w)] += 1
    ok(cnt[1:11] == row, f"table n={n}: {cnt[1:11]}")
# source Lemma 3 bound IB_k(n,t)/k^n <= k^(1-t/2) fails at k=2, n=11, t=9
ok(Fraction(192, 2**11) ** 2 > Fraction(1, 2**7), "lemma-3 counterexample")

# --- unbordered words -----------------------------------------------------------------
def unbordered(k, J):
    u = [0] * (J + 1)
    for j in range(1, J + 1):
        u[j] = k if j == 1 else (k * u[j - 1] - (u[j // 2] if j % 2 == 0 else 0))
    return u
def bordered(w):
    return any(w[:l] == w[len(w) - l:] for l in range(1, len(w)))
for k, J in ((2, 14), (3, 9)):
    u = unbordered(k, J)
    for j in range(1, J + 1):
        ok(u[j] == sum(1 for w in itertools.product(range(k), repeat=j) if not bordered(w)), f"u k={k} j={j}")

# --- renewal tail Pr(R_n >= r) = sum over j_1..j_r, 2*sum <= n of prod q_j (exact) ------
def tail_formula(k, n, r):
    u = unbordered(k, n)
    q = [Fraction(0)] + [Fraction(u[j], k ** (2 * j)) for j in range(1, n + 1)]
    # f[s] = sum over r-tuples with sum s of prod q
    f = {0: Fraction(1)}
    for _ in range(r):
        g = {}
        for s, v in f.items():
            for j in range(1, n + 1):
                if 2 * (s + j) <= n:
                    g[s + j] = g.get(s + j, 0) + v * q[j]
        f = g
    return sum(f.values(), Fraction(0))
for k, n in ((2, 12), (3, 8)):
    words = list(itertools.product(range(k), repeat=n))
    for r in range(1, n // 2 + 1):
        emp = Fraction(sum(1 for w in words if maxw_greedy(w) >= 2 * r), len(words))
        # R_n >= r  <=>  width >= 2r
        ok(emp == tail_formula(k, n, r), f"tail k={k} n={n} r={r}")

# --- A274190 triangle and A274199 -------------------------------------------------------
def triangle(Nmax):
    g = [[0] * (2 * Nmax + 2) for _ in range(Nmax + 1)]
    for n in range(Nmax + 1):
        g[n][0] = 1
        if n > 0:
            for m in range(1, n + 1):
                g[n][m] = g[n - 1][m - 1] + g[n - 1][2 * m]
    return g
g = triangle(160)
b190 = [int(l.split()[1]) for l in open(sys.argv[1]) if l.strip() and not l.startswith('#')]
flat = [g[n][m] for n in range(160) for m in range(n + 1)][:len(b190)]
ok(flat == b190[:len(flat)], "A274190 b-file")
b199 = [int(l.split()[1]) for l in open(sys.argv[2]) if l.strip() and not l.startswith('#')]
for d in range(0, 75):
    vals = {g[N][N - d] for N in range(2 * d, 160)}
    ok(len(vals) == 1, f"stabilization d={d}")
    if d < len(b199):
        ok(g[2 * d][d] == b199[d], f"A274199 d={d}")
ok(g[2 * 9 - 1][9 - 9] == 1, "sanity")

# compositions with a_{i+1} < 2 a_i
def comp_count(D):
    # c[s][last] = number of compositions of s with last part `last`
    c = [[0] * (D + 1) for _ in range(D + 1)]
    for s in range(1, D + 1):
        c[s][s] = 1
        for last in range(1, s):
            c[s][last] = sum(c[s - last][p] for p in range(1, D + 1) if last < 2 * p)
    return [1] + [sum(c[s]) for s in range(1, D + 1)]
A = comp_count(60)
for d in range(0, 61):
    ok(A[d] == g[2 * d][d], f"compositions d={d}")

# --- signed doubling chains: A(x) = 1/(1 - Q(x)) as formal series ------------------------
D = 60
p = [None] * (D + 1)  # p[j] polynomial (list of coefficients up to degree D)
def padd(a, b, s=1):
    return [a[i] + s * b[i] for i in range(D + 1)]
def pshift(a, i):
    return [0] * i + a[:D + 1 - i]
for j in range(1, D + 1):
    acc = [1] + [0] * D
    for i in range(1, j // 2 + 1):
        acc = padd(acc, pshift(p[i], i), -1)
    p[j] = acc
Q = [0] * (D + 1)
for j in range(1, D + 1):
    Q = padd(Q, pshift(p[j], j))
# check p_j(1/k) = u_j(k)/k^j for k = 2, 3
for k in (2, 3):
    u = unbordered(k, D)
    for j in range(1, 25):
        ok(sum(Fraction(c, k ** e) for e, c in enumerate(p[j])) == Fraction(u[j], k ** j), f"p_j k={k} j={j}")
inv = [1] + [0] * D  # 1/(1-Q)
for n in range(1, D + 1):
    inv[n] = sum(Q[i] * inv[n - i] for i in range(1, n + 1))
ok(inv == A[:D + 1], "A = 1/(1-Q)")

# --- the conjectured value versus 2/(1-q) - 1 ---------------------------------------------
getcontext().prec = 60
Abig = comp_count(400)
for k in range(2, 9):
    u = unbordered(k, 400)
    q = sum(Decimal(u[j]) / Decimal(k) ** (2 * j) for j in range(1, 401))
    lhs = 2 / (1 - q) - 1
    rhs = 1 + 2 * sum(Decimal(Abig[d]) / Decimal(k) ** d for d in range(1, 401))
    ok(abs(lhs - rhs) < Decimal(10) ** -30, f"E_k k={k}: {lhs} vs {rhs}")
    if k == 2:
        print("E_2 =", str(lhs)[:22])
    # exact E_{n,k} for moderate n via the renewal law (odd n: no empty centre)

# source Table 3 (four decimals; the table mixes rounding and truncation, e.g. E_8 = 1.32726... is
# printed 1.3272 while E_9 = 1.28168... is printed 1.2817, so allow one unit in the last place)
src_tab = {2: "6.4686", 3: "2.5908", 4: "1.9080", 5: "1.6314", 6: "1.4827", 7: "1.3902",
           8: "1.3272", 9: "1.2817", 10: "1.2472", 100: "1.0204"}
for k, s in src_tab.items():
    u = unbordered(k, 400)
    q = sum(Decimal(u[j]) / Decimal(k) ** (2 * j) for j in range(1, 401))
    ok(abs((2 / (1 - q) - 1) - Decimal(s)) <= Decimal("0.0001"), f"table3 k={k}")

# exact E_{n,2} by enumeration for n <= 16 versus the IB recurrence route (sanity on monotone approach)
for n in (15, 16):
    tot = sum(maxw_greedy(w) for w in itertools.product(range(2), repeat=n))
    print(f"E_{{{n},2}} =", float(Fraction(tot, 2 ** n)))

print("bad=", bad)
sys.exit(1 if bad else 0)
