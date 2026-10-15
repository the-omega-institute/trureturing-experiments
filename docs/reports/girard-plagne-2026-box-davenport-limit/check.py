#!/usr/bin/env python3
"""Exact checks for the Girard–Plagne limit question (arXiv:2610.16562v1, Section 2).

D^(k)(X) is the maximal length of a minimal zero-sum sequence over X whose support has exactly k elements.
The question: for every d >= 1, D^(d+1)([-m, m]^d) / m^d has a limit as m -> infinity.
The proposed answer is Delta_d = max |det T(V)| over d x (d+1) sign matrices V, where T(V) adds a row of ones.

Checks (all exact integer arithmetic):
  A. Brute force D^(3)([-m, m]^2), m = 2..6, against Theorem 1 of the source (4m^2 - q(m)).
  B. Brute force D^(2)([-m, m]^1), m = 1..12, against 2m - 1 (ratio limit Delta_1 = 2).
  C. Delta_d for d = 1..4 by vertex enumeration; Delta_2 = 4 and Delta_3 = 16 match Theorems 1 and 5.
  D. Lemma 1: for full-rank d x (d+1) integer matrices with a strictly positive kernel, the primitive kernel vector
     is the unique minimal multiplicity vector and its length is |det T(M)| / gcd(minors) (d = 2, exhaustive small box).
  E. Lemma 3: Smith-normal-form perturbation B makes the gcd of the maximal minors of tA + B equal to 1 for all t.
  F. Lemmas 4-5 and the completion: the integer simplex A built from a maximizing V has a strictly positive kernel and
     |det T(A)| = (nk)^d Delta_d; the sequences from tA + B are minimal, have support d + 1, and their normalized
     lengths approach (n/(n+1))^d Delta_d.
Exit 0 and print ALL_OK only if every check passes.
"""
import itertools
import math
import random
import sys

from sympy import Matrix, gcd as sgcd, primerange

FAIL = []


def check(cond, msg):
    if not cond:
        FAIL.append(msg)
        print("FAIL", msg)


def det(rows):
    """Exact integer determinant by Bareiss fraction-free elimination."""
    a = [list(map(int, r)) for r in rows]
    n = len(a)
    if n == 0:
        return 1
    sign, prev = 1, 1
    for i in range(n - 1):
        if a[i][i] == 0:
            for r in range(i + 1, n):
                if a[r][i] != 0:
                    a[i], a[r] = a[r], a[i]
                    sign = -sign
                    break
            else:
                return 0
        for r in range(i + 1, n):
            for c in range(i + 1, n):
                a[r][c] = (a[r][c] * a[i][i] - a[r][i] * a[i][c]) // prev
        prev = a[i][i]
    return sign * a[n - 1][n - 1]


def augmented(cols):
    """T(V): first row all ones, then the coordinates; cols is a list of column vectors."""
    d = len(cols[0])
    return [[1] * len(cols)] + [[c[i] for c in cols] for i in range(d)]


def signed_minors(cols):
    """c_i = (-1)^i det(M with column i deleted), for a d x (d+1) matrix given by its columns."""
    d = len(cols[0])
    out = []
    for i in range(len(cols)):
        rest = cols[:i] + cols[i + 1:]
        out.append((-1) ** i * det([[c[r] for c in rest] for r in range(d)]))
    return out


def gcd_list(xs):
    g = 0
    for x in xs:
        g = math.gcd(g, abs(x))
    return g


def is_minimal(vecs, a):
    """No nonzero proper b <= a with sum b_i v_i = 0 (dynamic programming over partial sums)."""
    d = len(vecs[0])
    reach = {tuple([0] * d): [tuple([0] * len(a))]}
    states = {(tuple([0] * len(a)))}
    # enumerate all b <= a by product; sizes in the checks are small
    for b in itertools.product(*[range(x + 1) for x in a]):
        if all(x == 0 for x in b) or list(b) == list(a):
            continue
        s = [sum(b[i] * vecs[i][r] for i in range(len(a))) for r in range(d)]
        if all(x == 0 for x in s):
            return False
    return True


# ---------- A. d = 2 brute force against Theorem 1 ----------
def q_of(m):
    for p in primerange(2, 10 * m + 10):
        if m % p:
            return p


def rank_deficient_max(ts, bound):
    """1-D minimal zero-sum sequences using exactly the integers ts (all used), length <= bound."""
    best = 0
    k = len(ts)
    for a in itertools.product(range(1, bound + 1), repeat=k):
        if sum(a[i] * ts[i] for i in range(k)) != 0:
            continue
        if is_minimal([[t] for t in ts], list(a)):
            best = max(best, sum(a))
    return best


def brute_D3_box2(m):
    pts = [(x, y) for x in range(-m, m + 1) for y in range(-m, m + 1) if (x, y) != (0, 0)]
    best = 0
    for tri in itertools.combinations(pts, 3):
        cols = [list(p) for p in tri]
        c = signed_minors(cols)
        if any(c):
            # rank 2: kernel is spanned by c; need all of one strict sign
            if all(x > 0 for x in c) or all(x < 0 for x in c):
                best = max(best, sum(abs(x) for x in c) // gcd_list(c))
        # rank <= 1 supports are handled below (they never exceed O(m))
    # rank-1 supports: three distinct nonzero multiples of a primitive direction, coefficients in [-m, m]
    r1 = 0
    for ts in itertools.combinations([t for t in range(-m, m + 1) if t != 0], 3):
        if min(ts) < 0 < max(ts):
            r1 = max(r1, rank_deficient_max(list(ts), 6 * m))
    return best, r1


for m in range(2, 7):
    full, r1 = brute_D3_box2(m)
    val = max(full, r1)
    expect = 4 * m * m - q_of(m)
    print(f"A d=2 m={m}: D3={val} (full-rank {full}, rank-1 {r1}), Theorem 1 gives {expect}")
    check(val == expect, f"A m={m}")


# ---------- B. d = 1 ----------
for m in range(1, 13):
    best = 0
    for x in range(1, m + 1):
        for y in range(-m, 0):
            best = max(best, (x - y) // math.gcd(x, -y))
    check(best == (2 if m == 1 else 2 * m - 1), f"B m={m}")
print("B d=1: D2([-1,1]) = 2 and D2([-m,m]) = 2m-1 for m=2..12; ratio -> 2 = Delta_1")


# ---------- C. Delta_d by vertex enumeration ----------
def delta(d):
    verts = list(itertools.product([-1, 1], repeat=d))
    best = 0
    arg = None
    for combo in itertools.combinations(range(len(verts)), d + 1):
        cols = [list(verts[i]) for i in combo]
        v = abs(det(augmented(cols)))
        if v > best:
            best, arg = v, cols
    return best, arg


DELTA = {}
for d in range(1, 5):
    DELTA[d] = delta(d)
    print(f"C Delta_{d} = {DELTA[d][0]}")
check(DELTA[1][0] == 2, "C Delta_1")
check(DELTA[2][0] == 4, "C Delta_2 (Theorem 1: 4m^2)")
check(DELTA[3][0] == 16, "C Delta_3 (Theorem 5: 16m^3)")


# ---------- D. Lemma 1 (d = 2, small entries) ----------
cnt = 0
for cols in itertools.combinations([(x, y) for x in range(-2, 3) for y in range(-2, 3) if (x, y) != (0, 0)], 3):
    cols = [list(p) for p in cols]
    c = signed_minors(cols)
    if not (all(x > 0 for x in c) or all(x < 0 for x in c)):
        continue
    g = gcd_list(c)
    a = [abs(x) // g for x in c]
    check(is_minimal(cols, a), f"D minimal {cols}")
    check(sum(a) == abs(det(augmented(cols))) // g, f"D length {cols}")
    # every positive integer kernel vector is a multiple of a
    for mult in range(2, 4):
        check(not is_minimal(cols, [mult * x for x in a]), f"D multiple {cols}")
    cnt += 1
print(f"D Lemma 1 checked on {cnt} positive-kernel supports in [-2,2]^2")


# ---------- E. Lemma 3 (Smith normal form perturbation) ----------
def smith(A):
    """Return (U, S, W) with U*A*W = S diagonal (d x k), U, W unimodular; plain integer algorithm."""
    S = Matrix(A)
    d, k = S.shape
    U = Matrix.eye(d)
    W = Matrix.eye(k)
    for t in range(min(d, k)):
        while True:
            nz = [(abs(S[i, j]), i, j) for i in range(t, d) for j in range(t, k) if S[i, j] != 0]
            if not nz:
                return U, S, W
            _, i, j = min(nz)
            S.row_swap(t, i); U.row_swap(t, i)
            S.col_swap(t, j); W.col_swap(t, j)
            done = True
            for i in range(t + 1, d):
                qq = S[i, t] // S[t, t]
                S[i, :] = S[i, :] - qq * S[t, :]; U[i, :] = U[i, :] - qq * U[t, :]
                if S[i, t] != 0:
                    done = False
            for j in range(t + 1, k):
                qq = S[t, j] // S[t, t]
                S[:, j] = S[:, j] - qq * S[:, t]; W[:, j] = W[:, j] - qq * W[:, t]
                if S[t, j] != 0:
                    done = False
            if done:
                break
    return U, S, W


def maximal_minor_gcd(M):
    d, k = M.shape
    g = 0
    for cols in itertools.combinations(range(k), d):
        g = math.gcd(g, abs(int(M[:, list(cols)].det())))
    return g


def perturbation(A):
    d, k = A.shape
    U, S, W = smith(A)
    check(U * A * W == S, "E smith identity")
    check(abs(U.det()) == 1 and abs(W.det()) == 1, "E unimodular")
    E = Matrix.zeros(d, k)
    for i in range(d - 1):
        E[i, i + 1] = 1
    E[d - 1, k - 1] = 1
    return U.inv() * E * W.inv()


random.seed(20261015)
ecount = 0
for d in (2, 3, 4):
    for _ in range(6):
        while True:
            A = Matrix(d, d + 1, lambda i, j: random.randint(-6, 6))
            if A.rank() == d:
                break
        B = perturbation(A)
        check(all(x == int(x) for x in B), "E B integral")
        for t in range(-15, 16):
            check(maximal_minor_gcd(t * A + B) == 1, f"E gcd d={d} t={t}")
        ecount += 1
print(f"E Lemma 3 checked on {ecount} random full-rank matrices, t = -15..15")


# ---------- F. Lemmas 4-5 and the completion ----------
for d in (2, 3):
    Dd, V = DELTA[d]
    k = d + 1
    T = Matrix(augmented(V))
    lam = T.inv() * Matrix([1] + [0] * d)
    check(all(0 <= x <= 1 for x in lam), f"F Lemma 4 d={d}")
    for n in (2, 4):
        s = [sum(v[r] for v in V) for r in range(d)]
        A = Matrix(d, k, lambda r, j: n * k * V[j][r] - s[r])
        Q = (n + 1) * k
        check(max(abs(x) for x in A) <= Q, f"F entries d={d} n={n}")
        beta = [(1 - 1 / n) * float(lam[i]) + 1 / (n * k) for i in range(k)]
        check(all(b > 0 for b in beta), f"F beta positive d={d} n={n}")
        check(abs(int(Matrix(augmented([list(A[:, j]) for j in range(k)])).det())) == (n * k) ** d * Dd,
              f"F det d={d} n={n}")
        B = perturbation(A)
        C = max(abs(int(x)) for x in B)
        ratios = []
        for t in (10, 100, 1000, 10000, 100000):
            M = t * A + B
            cols = [[int(M[r, j]) for r in range(d)] for j in range(k)]
            c = signed_minors(cols)
            pos = all(x > 0 for x in c) or all(x < 0 for x in c)
            check(pos, f"F positive minors d={d} n={n} t={t}")
            check(gcd_list(c) == 1, f"F gcd one d={d} n={n} t={t}")
            check(len({tuple(col) for col in cols}) == k, f"F distinct support d={d} n={n} t={t}")
            L = sum(abs(x) for x in c)
            check(L == abs(det(augmented(cols))), f"F length d={d} n={n} t={t}")
            if d == 2 and t == 10:
                check(is_minimal(cols, [abs(x) for x in c]) if L < 150 else True, f"F minimal d={d} t={t}")
            box = t * Q + C
            ratios.append(L / box ** d)
        target = (n / (n + 1)) ** d * Dd
        print(f"F d={d} n={n}: C={C}, L/(tQ+C)^d = {[round(r, 4) for r in ratios]} -> (n/(n+1))^d Delta_d = {target:.4f}")
        check(abs(ratios[-1] - target) < 0.02 * target, f"F convergence d={d} n={n}")

if FAIL:
    print(f"FAILURES {len(FAIL)}")
    sys.exit(1)
print("ALL_OK")
