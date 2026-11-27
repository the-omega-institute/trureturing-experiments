"""Exact minus-sector scales and a numerical rank scan for general cutoff-56 completions.

Exact part.  H is the identity except on the three isolated near edges, where it is
[[1, 1], [1, -1]].  S -> H^T S H maps the 37-dimensional space of free corrections onto
itself, so W = Delta - S and H^T W H correspond one to one with rank and positivity
preserved.  In this basis the three minus vectors m = v_1 - v_2 have fixed cross inner
products sigma, -sigma, sigma with one common sigma, fixed inner products of common
absolute value tau with the other plus vectors, fixed inner products below 7e-11 with
the other thirteen points, while the largest off-diagonal entry of Delta lies in
(3e-7, 4e-7).

Numerical part (numerical exploration only, real symmetric completions, two fixed seeds
per rank, no certificate).  Levenberg-Marquardt on the sum of squares of the smallest
19 - r eigenvalues of W = Delta - S over the 37 free entries.
"""
import sys
from itertools import combinations
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import fmpq, fmpq_mat

sys.set_int_max_str_digits(0)
prev = run_path(str(Path(__file__).with_name("ch39-symmetric-gram-lower-bound.py")))
PAIRS = ((0, 2), (1, 3), (17, 18))          # plus vector stored at p, minus vector at q
LARGE = (6, 7, 16, 8, 9, 14, 11)
MINUS = tuple(q for _, q in PAIRS)
PLUS = tuple(p for p, _ in PAIRS)


def exact_pm(data):
    delta = data["delta"]
    h = fmpq_mat(19, 19)
    for i in range(19):
        h[i, i] = 1
    for p, q in PAIRS:
        h[p, p], h[q, p], h[p, q], h[q, q] = 1, 1, 1, -1
    return h.transpose()*delta*h


def free_positions(dist):
    pair = {p: q for p, q in PAIRS} | {q: p for p, q in PAIRS}
    def free(i, j):
        if i == j or pair.get(i) == j:
            return True
        if i in pair or j in pair:
            return False
        return dist[i][j] <= 56
    return [(i, j) for i in range(19) for j in range(i, 19) if free(i, j)], free


def basis_change_is_bijective(positions):
    """H^T E H stays in the span of the free generators E, and the 37 images are independent."""
    h = fmpq_mat(19, 19)
    for i in range(19):
        h[i, i] = 1
    for p, q in PAIRS:
        h[p, p], h[q, p], h[p, q], h[q, q] = 1, 1, 1, -1
    allowed = set(positions) | {(j, i) for i, j in positions}
    rows = []
    for i, j in positions:
        e = fmpq_mat(19, 19)
        e[i, j] = e[j, i] = 1
        g = h.transpose()*e*h
        assert all(g[u, v] == 0 for u in range(19) for v in range(19) if (u, v) not in allowed)
        rows.append([g[u, v] for u, v in positions])
    assert fmpq_mat(rows).rank() == len(positions) == 37


def exact_facts(data, dpm, free):
    m = MINUS
    sig = dpm[m[0], m[1]]
    assert sig > 0 and dpm[m[1], m[2]] == sig and dpm[m[0], m[2]] == -sig
    assert fmpq(1320, 10**17) < sig < fmpq(1321, 10**17)
    tau = abs(dpm[2, 1])
    assert all(abs(dpm[q, p]) == tau for q in m for p in PLUS if not free(min(p, q), max(p, q)))
    assert fmpq(1322, 10**17) < tau < fmpq(1323, 10**17)
    others = [j for j in range(19) if j not in PLUS+MINUS]
    big = max(abs(dpm[q, j]) for q in m for j in others)
    assert fmpq(6, 10**11) < big < fmpq(7, 10**11)
    gram = max(abs(data["delta"][i, j]) for i, j in combinations(range(19), 2))
    assert fmpq(3, 10**7) < gram < fmpq(4, 10**7)
    return sig, tau, big, gram


class Scan:
    def __init__(self, dpm, free):
        s = np.full(19, 2048.0)
        s[list(LARGE)] = 8192.0
        for q in MINUS:
            s[q] = 2048.0/(2*1.95e-4)
        a = np.array([[float(dpm[i, j]) for j in range(19)] for i in range(19)])
        self.A = a*np.outer(s, s)
        self.free = free
        b = np.zeros((len(free), 19, 19))
        for k, (i, j) in enumerate(free):
            b[k, i, j] = b[k, j, i] = 1
        self.B = b

    def W(self, x):
        return self.A-np.tensordot(x, self.B, 1)

    def obj(self, x, m):
        return float((np.linalg.eigvalsh(self.W(x))[:m]**2).sum())

    def lm(self, x, m, iters, lam=1e-2):
        c = self.obj(x, m)
        iu = np.triu_indices(m)
        for _ in range(iters):
            w = self.W(x)
            _, v = np.linalg.eigh(w)
            u = v[:, :m]
            r = (u.T@w@u)[iu]
            jm = -np.einsum("ia,kij,jb->kab", u, self.B, u)[:, iu[0], iu[1]].T
            while True:
                dx = np.linalg.solve(jm.T@jm+lam*np.eye(len(x)), -jm.T@r)
                cn = self.obj(x+dx, m)
                if cn < c:
                    x, c, lam = x+dx, cn, max(lam/3, 1e-16)
                    break
                lam *= 5
                if lam > 1e12:
                    return x, c
            if c < 1e-28:
                break
        return x, c

    def start(self, seed):
        rng = np.random.default_rng(seed)
        x = np.zeros(len(self.free))
        for k, (i, j) in enumerate(self.free):
            if i == j:
                x[k] = self.A[i, i]-rng.uniform(0.5, 2.0) if i in MINUS else self.A[i, i]*rng.uniform(0, 0.9)
            else:
                x[k] = self.A[i, j]-rng.normal()*0.3 if i in PLUS+MINUS else rng.normal()*0.3
        return x


def verify():
    data = prev["setup"]()
    dpm = exact_pm(data)
    free, isfree = free_positions(data["dist"])
    assert len(free) == 37
    basis_change_is_bijective(free)
    sig, tau, big, gram = exact_facts(data, dpm, isfree)
    print("minus-minus sigma", float(sig), "minus-plus", float(tau), "minus-other <=", float(big),
          "max old Gram entry", float(gram))
    scan = Scan(dpm, free)
    # Rank 12, two fixed seeds.
    for seed in (0, 1):
        _, c = scan.lm(scan.start(seed), 7, 3000)
        assert c < 1e-24
        print("rank 12 seed", seed, "objective", c)
    # Rank 10, two fixed seeds.
    for seed in (0, 1):
        _, c = scan.lm(scan.start(seed), 9, 1500)
        assert c > 0.3
        print("rank 10 seed", seed, "objective", c)
    # Rank 11, two fixed seeds: objective after 3000 and 6000 steps, minus norms.
    grown = []
    for seed in (0, 1):
        x, c1 = scan.lm(scan.start(seed), 8, 3000)
        w1 = scan.W(x)
        x, c2 = scan.lm(x, 8, 3000)
        w2 = scan.W(x)
        assert 5e-7 < c2 < c1 < 5e-6 and 1.1e-6 < c2 < 1.2e-6
        ratio = {q: w2[q, q]/w1[q, q] for q in MINUS}
        up = {q for q in MINUS if ratio[q] > 1.2}
        assert len(up) == 2 and all(ratio[q] < 1.3 for q in up)
        assert all(abs(ratio[q]-1) < 0.01 for q in MINUS if q not in up)
        grown.append(up)
        print("rank 11 seed", seed, "objective", c1, "->", c2, "minus norms",
              [round(w1[q, q], 1) for q in MINUS], "->", [round(w2[q, q], 1) for q in MINUS])
    assert grown[0] != grown[1]


if __name__ == "__main__":
    verify()
