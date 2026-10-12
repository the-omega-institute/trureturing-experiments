"""Numerical check of Hanaoka-Konno Conjecture 4.2 (arXiv:2010.04518v3, QIC 21 (2021) 409-422) for the Riesz walk.

The walk is built literally from the paper's Section 2: Riesz moments -> Caratheodory function F -> Schur function f
-> Schur algorithm (Verblunsky parameters alpha_k) -> the block matrix U^(s) -> Psi_{t+1} = U^(s) Psi_t with
Psi_0 = [1,0]^T delta_0, and mu_t(x) = |Psi_t^L(x)|^2 + |Psi_t^R(x)|^2.
Conjecture 4.2 (equivalent form, paper p.12): mu_{2t}(0) = mu_{8t}(0) and
mu_{2t}(2x-1) + mu_{2t}(2x) = sum_{y=8x-7}^{8x} mu_{8t}(y) for t, x >= 1.
"""
import mpmath as mp
mp.mp.dps = 60
TMAX = 32                 # check t = 1..TMAX, so times up to 8*TMAX
NT = 8 * TMAX
NA = 2 * NT + 8           # Verblunsky parameters needed
DEG = NA + 8              # power-series degree

# Riesz moments: expand prod_{k>=1} (1 + (z^{4^k} + z^{-4^k})/2) far enough.
mom = {0: mp.mpf(1)}
k = 1
while 4 ** k <= 4 * DEG:
    a = 4 ** k; new = {}
    for e, c in mom.items():
        for d, w in ((0, 1), (a, mp.mpf(1)/2), (-a, mp.mpf(1)/2)):
            new[e + d] = new.get(e + d, 0) + c * w
    mom = {e: c for e, c in new.items() if abs(e) <= 4 * DEG}
    k += 1
mu = lambda j: mom.get(j, mp.mpf(0))

# Power series helpers (lists of coefficients, length DEG).
def mul(a, b):
    r = [mp.mpf(0)] * DEG
    for i, x in enumerate(a):
        if x == 0: continue
        for j in range(DEG - i):
            if b[j] != 0: r[i + j] += x * b[j]
    return r
def inv(a):
    r = [mp.mpf(0)] * DEG; r[0] = 1 / a[0]
    for n in range(1, DEG):
        s = mp.mpf(0)
        for i in range(1, n + 1):
            if a[i] != 0: s += a[i] * r[n - i]
        r[n] = -s / a[0]
    return r

F = [mp.mpf(1)] + [2 * mu(n) for n in range(1, DEG)]       # moments are real
num = [F[0] - 1] + F[1:]; den = [F[0] + 1] + F[1:]
q = mul(num, inv(den))                                     # (F-1)/(F+1), q[0] = 0
f = q[1:] + [mp.mpf(0)]                                    # f = q / z
alpha = []
fk = f[:]
for kk in range(NA):
    a0 = fk[0]; alpha.append(a0)
    numk = [fk[0] - a0] + fk[1:]
    denk = [1 - a0 * fk[0]] + [-a0 * c for c in fk[1:]]
    g = mul(numk, inv(denk))
    fk = g[1:] + [mp.mpf(0)]
nz = [i for i, a in enumerate(alpha) if abs(a) > mp.mpf(10) ** -40]
assert all(i % 4 == 3 for i in nz), nz[:10]
print("nonzero alpha indices (first 6):", nz[:6], "values:", [mp.nstr(alpha[i], 8) for i in nz[:4]])
rho = [mp.sqrt(1 - a * a) for a in alpha]

# U^(s) on positions 0..P-1, coordinates (x, L/R) -> 2x, 2x+1.
P = NT + 2
def blocks():
    A, R = alpha, rho
    Rb = {0: [[A[0], R[0]], [R[0] * A[1], -A[0] * A[1]]]}
    for x in range(1, P):
        Rb[x] = [[-A[2*x-1] * A[2*x], -A[2*x-1] * R[2*x]], [R[2*x] * A[2*x+1], -A[2*x] * A[2*x+1]]]
    Pb = {x: [[R[2*x] * R[2*x+1], -A[2*x] * R[2*x+1]], [0, 0]] for x in range(P)}
    Qb = {x: [[0, 0], [R[2*x-1] * A[2*x], R[2*x-1] * R[2*x]]] for x in range(1, P)}
    return Rb, Pb, Qb
Rb, Pb, Qb = blocks()
def step(psi):
    new = [mp.mpf(0)] * (2 * P)
    for x in range(P):
        for (blk, src) in ((Rb.get(x), x), (Pb.get(x - 1), x - 1), (Qb.get(x + 1), x + 1)):
            if blk is None or src < 0 or src >= P: continue
            for i in range(2):
                for j in range(2):
                    if blk[i][j] != 0: new[2*x+i] += blk[i][j] * psi[2*src+j]
    return new
psi = [mp.mpf(0)] * (2 * P); psi[0] = mp.mpf(1)
dist = {}
for n in range(NT + 1):
    dist[n] = [psi[2*x] ** 2 + psi[2*x+1] ** 2 for x in range(P)]
    psi = step(psi)
tot_err = max(abs(sum(dist[n]) - 1) for n in dist)
print("max |total probability - 1| =", mp.nstr(tot_err, 3))
print("return probabilities mu_4, mu_16, mu_64, mu_256:", [mp.nstr(dist[n][0], 10) for n in (4, 16, 64, 256) if n <= NT])

bad = 0; worst = mp.mpf(0)
for t in range(1, TMAX + 1):
    d2, d8 = dist[2 * t], dist[8 * t]
    e = abs(d2[0] - d8[0]); worst = max(worst, e)
    if e > mp.mpf(10) ** -30: bad += 1; print("FAIL origin t=", t, mp.nstr(e, 5))
    for x in range(1, t + 1):
        lhs = d2[2*x-1] + d2[2*x]
        rhs = sum(d8[y] for y in range(8*x - 7, 8*x + 1))
        e = abs(lhs - rhs); worst = max(worst, e)
        if e > mp.mpf(10) ** -30: bad += 1; print("FAIL t=", t, "x=", x, mp.nstr(e, 5))
print("checked t = 1 ..", TMAX, "worst discrepancy", mp.nstr(worst, 3))
# Control: the identity fails for the factor 2 (time 2t vs 4t) at some t, so it is not vacuous.
ctrl = 0
for t in range(1, TMAX // 2 + 1):
    d2, d4 = dist[2 * t], dist[4 * t]
    for x in range(1, t + 1):
        if abs((d2[2*x-1] + d2[2*x]) - sum(d4[y] for y in range(4*x - 3, 4*x + 1))) > 1e-6: ctrl += 1
print("control (factor 2 instead of 4) mismatches:", ctrl)
if ctrl == 0: print("FAIL control vacuous"); bad += 1
print("bad=", bad)
raise SystemExit(1 if bad else 0)
