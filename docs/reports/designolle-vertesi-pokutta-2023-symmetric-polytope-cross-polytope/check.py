#!/usr/bin/env python3
"""Designolle–Vértesi–Pokutta, arXiv:2310.20677v3 (Phys. Rev. A 109, 022205), Section V: "the symmetric polytope
is affinely equivalent to the cross-polytope in dimension ceil(m/2) [...] We conjecture, however, that it holds in
general."  Checks that the conjecture fails for N = 3 parties and m = 17 inputs.

Conventions: inputs are 0..m-1 with e_m = -e_0 (antiperiodic), so that the GHZ tensor cos(pi/m sum x) of Eq. (7)
is invariant. A local strategy a in {+-1}^m is the polynomial A(X) = sum a_x X^x in R[X]/(X^m + 1); the
deterministic tensor of (a, b, c) has negacyclic total-index structure D = A B C. P(D) = (d_0, (d_r - d_{m-r})/2 for
1 <= r < m/2) are the invariant coordinates of Eq. (27).

(A) The source's Reynolds operator, built literally from the generators g1..g4 of Eq. (25) by closing the group of
    signed permutations, equals P(D)/m^2 on the coordinates p~_{r,0,0}: on all strategies for m = 4, 5 and on 3000
    random strategies for m = 6, 7.
(B) The number of distinct symmetrised deterministic points matches the source's Table II for N = 3, m = 3..9.
(C) For odd m = 5..13 every symmetrised point lies in conv{+-V_j}, V_j = P(S_j^3), S_j = sign cos(pi j x/m), j odd:
    the projected polytope is that cross-polytope for those m (consistent with the source's m <= 12 computations).
(D) m = 17: V (9 x 9 integer matrix) has nonzero determinant; |A(zeta^j)| <= csc(pi/34) with exactly 34 maximizers
    for odd j < 17 (2 for j = 17), so each +-V_j is a uniquely exposed vertex; the strategy (a, b, c) gives the point
    p = P(ABC) with w.p = 8421 > 8417 = max |w.V_j|. Hence the projected polytope has more than 18 vertices and is
    not affinely equivalent to the 9-dimensional cross-polytope."""
import sys, itertools, math, random
import numpy as np
import sympy as sp

bad = 0
def report(label, ok):
    global bad
    print(label, ok)
    bad += 0 if ok else 1

def negacyclic(a, b):
    m = len(a); out = [0] * m
    for x in range(m):
        if a[x] == 0: continue
        for y in range(m):
            k = x + y
            if k < m: out[k] += a[x] * b[y]
            else: out[k - m] -= a[x] * b[y]
    return out

def P(d):
    m = len(d)
    return tuple([sp.Integer(d[0])] + [sp.Rational(d[r] - d[m - r], 2) for r in range(1, (m + 1) // 2)])

# ---------------- (A) literal Reynolds operator ----------------
def norm_index(v, m):
    q, r = divmod(v, m)
    return (-1) ** q, r

def group_elements(m, N=3):
    basis = list(itertools.product(range(m), repeat=N)); pos = {x: i for i, x in enumerate(basis)}; n = len(basis)
    def make(f):
        perm = np.zeros(n, dtype=np.int64); sign = np.zeros(n, dtype=np.int64)
        for i, x in enumerate(basis):
            s, y = f(x); perm[i] = pos[y]; sign[i] = s
        return (perm, sign)
    def g1(x): return 1, (x[1], x[0]) + x[2:]
    def g2(x): return 1, (x[-1],) + x[:-1]
    def g3(x):
        s1, y1 = norm_index(x[0] + 1, m); s2, y2 = norm_index(x[1] - 1, m)
        return s1 * s2, (y1, y2) + x[2:]
    def g4(x):
        s = (-1) ** N; ys = []
        for v in x:
            t, y = norm_index(m - v, m); s *= t; ys.append(y)
        return s, tuple(ys)
    gens = [make(g) for g in (g1, g2, g3, g4)]
    ident = (np.arange(n), np.ones(n, dtype=np.int64))
    key = lambda g: (g[0].tobytes(), g[1].tobytes())
    seen = {key(ident): ident}; frontier = [ident]
    while frontier:
        nxt = []
        for h in frontier:
            for g in gens:
                comp = (g[0][h[0]], h[1] * g[1][h[0]])
                k = key(comp)
                if k not in seen: seen[k] = comp; nxt.append(comp)
        frontier = nxt
    return list(seen.values()), basis, pos

def reynolds(d, G):
    out = np.zeros(len(d), dtype=object)
    for perm, sign in G:
        v = np.zeros(len(d), dtype=object)
        v[perm] = sign * d
        out = out + v
    return out / len(G)

random.seed(3)
for m, exhaustive in [(4, True), (5, True), (6, False), (7, False)]:
    G, basis, pos = group_elements(m)
    strategies = (itertools.product(itertools.product([1, -1], repeat=m), repeat=3) if exhaustive
                  else ([tuple(random.choice([1, -1]) for _ in range(m)) for _ in range(3)] for _ in range(3000)))
    ok = True; count = 0
    # the Reynolds operator is linear; compute Gamma on each basis vector's orbit-average coordinate functional once
    coords = [pos[(r, 0, 0)] for r in range(m)]
    # functional: Gamma(d)_{(r,0,0)} = (1/|G|) sum_g (g.d)_{(r,0,0)} = sum_x d_x * c_r[x]
    cfun = []
    for c in coords:
        acc = np.zeros(len(basis), dtype=np.int64)
        for perm, sign in G:
            # (g.d)_c = sign[x] d_x where perm[x] = c
            x = np.where(perm == c)[0][0]; acc[x] += sign[x]
        cfun.append(acc)
    for a, b, c in strategies:
        d = np.array([a[x] * b[y] * c[z] for (x, y, z) in basis], dtype=np.int64)
        gam = [sp.Rational(int(f @ d), len(G)) for f in cfun]
        D = negacyclic(negacyclic(list(a), list(b)), list(c))
        pred = [v / m**2 for v in P(D)]
        # Eq. (27): p~_{r,0,0} = -p~_{m-r,0,0} for r > m/2, and p~_{m/2,0,0} = 0 for even m
        mirror = all(gam[r] == -gam[m - r] for r in range(m // 2 + 1, m)) and (m % 2 == 1 or gam[m // 2] == 0)
        if gam[:len(pred)] != pred or not mirror: ok = False
        count += 1
    report(f'(A) m={m}: |G|={len(G)}, {count} strategies, Reynolds coordinates equal P(ABC)/m^2:', ok)

# ---------------- (B) counts against the source's Table II (N = 3) ----------------
def all_signs(m): return np.array(list(itertools.product([1, -1], repeat=m)), dtype=np.int64)
def shift_stack(E):
    """S[e, y, k] = coefficient k of X^y * E (negacyclic)."""
    n, m = E.shape; S = np.zeros((n, m, m), dtype=np.int64)
    for y in range(m):
        S[:, y, y:] = E[:, :m - y]; S[:, y, :y] = -E[:, m - y:]
    return S
def products_E(m, a_rows):
    B = all_signs(m)
    out = []
    for a in a_rows:
        Sa = shift_stack(a[None, :])[0]          # Sa[y, k]: X^y * A
        out.append(B @ Sa)                       # E = A*B for all b
    return np.unique(np.concatenate(out), axis=0)
def Pmat(m):
    h = (m + 1) // 2 if m % 2 else m // 2
    M = np.zeros((h, m), dtype=np.float64); M[0, 0] = 1
    for r in range(1, h): M[r, r] = 0.5; M[r, m - r] = -0.5
    return M
table = {3: 10, 4: 10, 5: 60, 6: 100, 7: 640, 8: 1540, 9: 10032}
for m, expected in table.items():
    E = products_E(m, all_signs(m)); C = all_signs(m)
    S = shift_stack(E)
    pts = set()
    for i in range(0, len(E), 4096):
        D = np.einsum('cy,eyk->eck', C, S[i:i + 4096]).reshape(-1, m)
        Q = np.round(2 * D @ Pmat(m).T).astype(np.int64)
        pts.update(map(tuple, np.unique(Q, axis=0)))
    report(f'(B) m={m}: {len(pts)} distinct symmetrised points, Table II gives {expected}:', len(pts) == expected)

# ---------------- (C) containment in conv(+-V_j) for odd m ----------------
def sign_cos(j, m): return np.array([1 if math.cos(math.pi * j * x / m) > 0 else -1 for x in range(m)], dtype=np.int64)
def Vmatrix(m):
    cols = []
    for j in range(1, m + 1, 2):
        s = list(sign_cos(j, m)); D = negacyclic(negacyclic(s, s), s)
        cols.append([int(2 * v) for v in P(D)])   # 2 * P(D) is integral
    return np.array(cols, dtype=np.float64).T     # columns 2*V_j
def necklace_reps(m):
    seen = set(); reps = []
    for a in itertools.product([1, -1], repeat=m):
        if a in seen: continue
        reps.append(a); cur = list(a)
        for _ in range(2 * m):
            cur = [-cur[-1]] + cur[:-1]; seen.add(tuple(cur)); seen.add(tuple(-v for v in cur))
    return np.array(reps, dtype=np.int64)
for m in (5, 7, 9, 11, 13):
    V2 = Vmatrix(m); d = V2.shape[0]
    Vinv = np.linalg.inv(V2)
    E = products_E(m, necklace_reps(m))
    best = 0.0
    for sigma in itertools.product([1, -1], repeat=d):
        f = np.array(sigma, dtype=np.float64) @ Vinv @ (2 * Pmat(m))      # functional on D: sigma . V^-1 P(D)
        W = np.zeros((m, m))
        for j in range(m):
            for y in range(m):
                k = j + y
                W[j, y] = f[k] if k < m else -f[k - m]
        val = np.abs(E @ W).sum(axis=1).max()
        best = max(best, val)
    report(f'(C) m={m}: max over symmetrised points of ||V^-1 q||_1 = {best:.12f} <= 1:', best <= 1 + 1e-9)

# ---------------- (D) m = 17 certificate ----------------
m = 17
V = sp.Matrix([[sp.Integer(v) for v in P(negacyclic(negacyclic(list(sign_cos(j, m)), list(sign_cos(j, m))), list(sign_cos(j, m))))]
               for j in range(1, 18, 2)]).T
V_seat = sp.Matrix([[145,145,145,145,145,145,145,145,289],[143,127,95,47,-17,-73,-113,-137,-289],[137,73,-47,-127,-143,-95,17,113,289],
                    [127,-17,-137,-113,47,143,95,-73,-289],[113,-95,-127,73,137,-47,-143,17,289],[95,-137,-17,143,-73,-113,127,47,-289],
                    [73,-143,113,17,-127,137,-47,-95,289],[47,-113,143,-137,95,-17,-73,127,-289],[17,-47,73,-95,113,-127,137,-143,289]])
report('(D) V matches the displayed matrix:', V == V_seat)
detV = V.det(); print('(D) det V =', detV)
report('(D) det V = 98400527135576109875200 (nonzero):', detV == 98400527135576109875200)
R = 1 / math.sin(math.pi / 34)
B17 = all_signs(17)
for j in range(1, 18, 2):
    z = np.exp(1j * math.pi * j * np.arange(17) / 17)
    vals = np.abs(B17 @ z); mx = vals.max(); cnt = int(np.sum(vals > mx - 1e-9))
    s = sign_cos(j, 17)
    if j < 17:
        ok = abs(mx - R) < 1e-9 and cnt == 34
        rot = {tuple(s)}; cur = list(s)
        for _ in range(33):
            cur = [-cur[-1]] + cur[:-1]; rot.add(tuple(cur))
        ok = ok and set(map(tuple, B17[vals > mx - 1e-9])) == rot
    else:
        ok = abs(mx - 17) < 1e-9 and cnt == 2
    lv = sum(float(V[r, (j - 1) // 2]) * (1 if r == 0 else 2 * math.cos(math.pi * j * r / 17)) for r in range(9))
    target = R**3 if j < 17 else 17.0**3
    report(f'(D) j={j}: max |A(zeta^j)| = {mx:.10f} with {cnt} maximizers (the signed rotations of s_j); l_j(V_j) = {lv:.6f} vs {target:.6f}:',
           ok and abs(lv - target) < 1e-6)
a = [1,-1,1,-1,-1,1,-1,1,-1,1,-1,1,-1,-1,1,-1,1]
b = [-1,1,-1,1,-1,1,-1,1,1,-1,1,-1,1,-1,1,1,-1]
c = [1,-1,1,1,-1,1,-1,1,-1,1,-1,-1,1,-1,1,-1,1]
dd = negacyclic(negacyclic(a, b), c)
report('(D) ABC coefficients as displayed:', dd == [-141,133,-109,69,-13,-43,91,-123,139,-139,123,-91,43,13,-69,109,-133])
p = sp.Matrix(P(dd)); print('(D) p =', list(p))
w = sp.Matrix([-15,27,-21,8,59,24,-5,-6,10])
wV = list((w.T * V)); wp = (w.T * p)[0]
print('(D) w.V_j =', wV, ' w.p =', wp)
report('(D) max |w.V_j| = 8417 < w.p = 8421:', max(abs(t) for t in wV) == 8417 and wp == 8421)
print('bad=', bad)
sys.exit(1 if bad else 0)
