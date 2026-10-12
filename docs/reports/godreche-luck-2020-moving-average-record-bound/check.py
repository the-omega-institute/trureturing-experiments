"""Check of the proposed proof of the bound suggested by Godreche-Luck, 'Records for the moving average of a time
series' (arXiv:1907.07598v2, J. Stat. Mech. (2020) 023201), Section 6.6: 'the data shown in figure 1 suggest that
(6.33) provides an absolute upper bound for Q_n', i.e. Q_n <= 3/(2n) for n >= 3, where X_1..X_n are iid with a
continuous density, Y_i = X_{i-1} + X_i (2 <= i <= n) and Q_n = Pr(Y_n > max(Y_2, ..., Y_{n-1})).

(1) Deterministic counting identity behind the proof: for a configuration with distinct cyclic sums
    C_i = X_i + X_{i+1} (indices mod n), among the 2n oriented paths (delete one cyclic edge, choose a direction)
    exactly 2 + [the two largest cyclic sums are adjacent] end in a record.  Checked exactly on random integer
    configurations with distinct cyclic sums, n = 3..9.
(2) Monte Carlo of n Q_n for several distributions (exponential, uniform, lognormal, Pareto with index 1/2,
    half-Cauchy), n = 3..12, against 3/2; and of 2 n Q_n against 2 + Pr(two largest cyclic sums adjacent).
(3) Exploration for general window p (Section 7: 'expected ... to provide an absolute upper bound for Q_n^(p)
    for all n >= p'): Monte Carlo of n Q_n^(p) against R_p = (2p-1)!/(2^(2p-2)(p-1)!^2) for p = 3, 4.  Reported,
    not asserted: for p <= n <= 2p-2 the universal values (Q_p = 1, Q_(p+1) = 1/2, ...) already exceed R_p/n.
"""
import math
import numpy as np
rng = np.random.default_rng(20191106)
bad = 0

def record_count_oriented(x):
    n = len(x); C = [x[i] + x[(i + 1) % n] for i in range(n)]
    cnt = 0
    for cut in range(n):                       # delete cyclic edge `cut` (between positions cut and cut+1)
        order = [(cut + 1 + j) % n for j in range(n)]
        for seq in (order, order[::-1]):
            sums = [x[seq[j - 1]] + x[seq[j]] for j in range(1, n)]
            if all(sums[-1] > s for s in sums[:-1]): cnt += 1
    srt = sorted(range(n), key=lambda i: C[i], reverse=True)
    adj = (srt[1] - srt[0]) % n in (1, n - 1)
    return cnt, 2 + (1 if adj else 0)

tested = 0
for n in range(3, 10):
    for _ in range(3000):
        x = list(rng.integers(0, 10**6, size=n))
        C = [x[i] + x[(i + 1) % n] for i in range(n)]
        if len(set(C)) < n: continue
        a, b = record_count_oriented(x); tested += 1
        if a != b: bad += 1
print("(1) identity checked on", tested, "configurations, bad so far", bad)

dists = {
    "exponential": lambda s: rng.exponential(size=s),
    "uniform": lambda s: rng.random(size=s),
    "lognormal": lambda s: rng.lognormal(size=s),
    "pareto1/2": lambda s: (1 - rng.random(size=s)) ** (-2.0),
    "halfcauchy": lambda s: np.abs(rng.standard_cauchy(size=s)),
}
N = 400000
for name, gen in dists.items():
    row = []
    for n in range(3, 13):
        X = gen((N, n)); Y = X[:, 1:] + X[:, :-1]
        Q = np.mean(Y[:, -1] > Y[:, :-1].max(axis=1)) if n > 2 else 1.0
        C = X + np.roll(X, -1, axis=1)
        idx = np.argsort(-C, axis=1)[:, :2]
        adj = np.mean(((idx[:, 1] - idx[:, 0]) % n == 1) | ((idx[:, 0] - idx[:, 1]) % n == 1))
        se = math.sqrt(Q * (1 - Q) / N)
        if n * Q > 1.5 + 5 * n * se: bad += 1
        if abs(2 * n * Q - (2 + adj)) > 10 * 2 * n * se: bad += 1
        row.append(f"{n}:{n*Q:.4f}")
    print("(2)", name, " ".join(row))

def Rp(p): return math.factorial(2 * p - 1) / (2 ** (2 * p - 2) * math.factorial(p - 1) ** 2)
for p in (3, 4):
    for name in ("exponential", "pareto1/2", "halfcauchy"):
        row = []
        for n in range(p, p + 10):
            X = dists[name]((N, n)); W = np.lib.stride_tricks.sliding_window_view(X, p, axis=1).sum(axis=2)
            Q = np.mean(W[:, -1] > W[:, :-1].max(axis=1)) if W.shape[1] > 1 else 1.0
            row.append(f"{n}:{n*Q:.3f}")
        print(f"(3) p={p} R_p={Rp(p):.4f} {name}", " ".join(row))
print("bad=", bad)
