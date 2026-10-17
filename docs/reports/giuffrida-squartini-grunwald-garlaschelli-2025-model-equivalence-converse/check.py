"""Giuffrida–Squartini–Grünwald–Garlaschelli arXiv:2307.05645v4, Section IV.C: the reverse of implication (71),
"D_KL^NML = o(n) implies D_KL(c) = o(n) for all c", fails.

Family: binary n x n matrices G, scalar statistic c(G) = sum_i (1 + 2^-i) g_{i1}; observation = checkerboard.
1. Direct evaluation for n = 2, 4 (mpmath, 50 digits) over all 2^(n^2) matrices, from the source's definitions:
   (16) microcanonical law, (20) canonical law, (21) expectation-matched parameter, (22) measure divergence,
   NML laws from the maximized likelihoods, (69) model divergence. The canonical maximum likelihood of each matrix
   is found by maximizing over theta (concave in theta; the supremum is a limit for the two extreme fibers).
   Results agree with the fiber reduction of the issue (Lemmas 1–4).
2. Fiber reduction (numpy, float64, log-sum-exp) for even n <= 14: W_n <= 2n, 0 <= D_NML <= ln W_n
   (so D_NML/n -> 0), and n ln 2 - ln 4 <= D_KL(c*) <= n ln 2.
3. Lemma 5 at n = 4: the matched parameter equals -16 ln q_4, q_4 the root of q P'(q) = 42 P(q),
   P(q) = (1+q^24)(1+q^20)(1+q^18)(1+q^17).
"""
import itertools, sys
from collections import Counter
import numpy as np
from mpmath import mp, mpf, exp, log, findroot, diff

mp.dps = 50


def stat_exact(col):
    return sum((1 + mpf(2) ** (-(i + 1))) * b for i, b in enumerate(col))


def checkerboard_col(n):  # g_{i1} = 1 iff i + 1 is even (indices from 1)
    return tuple(1 if (i + 2) % 2 == 0 else 0 for i in range(n))


def mp_theta(values, mult, target):
    """Solve <c>_theta = target, eq. (21), by bisection (mean is decreasing in theta)."""
    def mean(t):
        z = sum(m * exp(-t * v) for v, m in zip(values, mult))
        return sum(m * v * exp(-t * v) for v, m in zip(values, mult)) / z
    lo, hi = mpf(-1), mpf(1)
    while mean(lo) < target: lo *= 2
    while mean(hi) > target: hi *= 2
    for _ in range(300):
        mid = (lo + hi) / 2
        if mean(mid) > target: lo = mid
        else: hi = mid
    return (lo + hi) / 2


def direct(n):
    mats = list(itertools.product((0, 1), repeat=n * n))
    st = [stat_exact(tuple(m[i * n] for i in range(n))) for m in mats]
    cnt = Counter(st); vals = sorted(cnt); mult = [cnt[v] for v in vals]
    lo_v, hi_v = vals[0], vals[-1]
    lcan = {}
    for v in vals:
        if v in (lo_v, hi_v):
            lcan[v] = 1 / mpf(cnt[v])          # sup_theta P_can(G; theta), attained in the limit theta -> +-inf
        else:
            t = mp_theta(vals, mult, v)
            z = sum(m * exp(-t * u) for u, m in zip(vals, mult))
            lcan[v] = exp(-t * v) / z
    Lmic = [1 / mpf(cnt[c]) for c in st]; Lcan = [lcan[c] for c in st]
    Zm, Zc = sum(Lmic), sum(Lcan)
    dnml = sum((a / Zm) * log((a / Zm) / (b / Zc)) for a, b in zip(Lmic, Lcan))
    cstar = stat_exact(checkerboard_col(n))
    t = mp_theta(vals, mult, cstar)
    z = sum(m * exp(-t * u) for u, m in zip(vals, mult))
    pmic = 1 / mpf(cnt[cstar])
    dkl = sum(pmic * log(pmic / (exp(-t * c) / z)) for c in st if c == cstar)
    return dnml, dkl, t


def reduced_np(n):
    cols = np.array(list(itertools.product((0, 1), repeat=n)), dtype=float)
    wts = 1 + 2.0 ** (-np.arange(1, n + 1))
    e = np.sort(cols @ wts); N = len(e)
    logw = np.zeros(N)
    for k in range(1, N - 1):
        t = 0.0
        for _ in range(200):
            a = -t * e; a -= a.max(); p = np.exp(a); p /= p.sum()
            mean = p @ e; var = p @ (e - mean) ** 2
            step = (mean - e[k]) / var       # d mean / d theta = -var
            t += step
            if abs(step) < 1e-13 * (1 + abs(t)): break
        a = -t * e; m = a.max()
        logw[k] = a[k] - (m + np.log(np.exp(a - m).sum()))
    W = np.exp(logw).sum()
    dnml = np.log(W) - np.mean(np.log(N) + logw)
    cstar = float(np.array(checkerboard_col(n)) @ wts)
    k = int(np.argmin(np.abs(e - cstar)))
    return dnml, -logw[k], W


def main():
    fails = []
    for n in (2, 4):
        dn, dk, t = direct(n)
        rn, rk, W = reduced_np(n)
        print(f"n={n}: direct D_NML={float(dn):.12f} D_KL(c*)={float(dk):.12f} theta*={float(t):.10f}; "
              f"fiber reduction D_NML={rn:.12f} D_KL(c*)={rk:.12f}")
        if abs(float(dn) - rn) > 1e-9 or abs(float(dk) - rk) > 1e-9: fails.append(("direct-vs-reduced", n))
    for n in range(2, 15, 2):
        dn, dk, W = reduced_np(n)
        lo, hi = n * np.log(2) - np.log(4), n * np.log(2)
        ok = -1e-12 <= dn <= np.log(W) + 1e-12 and W <= 2 * n and lo - 1e-9 <= dk <= hi + 1e-9
        print(f"n={n:2d}: W_n={W:.6f} (2n={2*n}) D_NML={dn:.6f} D_NML/n={dn/n:.6f} D_KL(c*)={dk:.6f} "
              f"D_KL/n={dk/n:.6f} bounds [{lo:.4f}, {hi:.4f}] {'OK' if ok else 'FAIL'}")
        if not ok: fails.append(("bounds", n))
    P = lambda q: (1 + q ** 24) * (1 + q ** 20) * (1 + q ** 18) * (1 + q ** 17)
    q4 = findroot(lambda q: q * diff(P, q) - 42 * P(q), mpf('1.01'))
    _, _, t4 = direct(4)
    print(f"Lemma 5 (n=4): q_4 = {float(q4):.15f}; -16 ln q_4 = {float(-16 * log(q4)):.15f}; matched theta = {float(t4):.15f}")
    if abs(-16 * log(q4) - t4) > mpf(10) ** -25: fails.append("lemma5")
    print("RESULT", "PASS" if not fails else f"FAIL {fails}")
    sys.exit(0 if not fails else 1)


if __name__ == "__main__":
    main()
