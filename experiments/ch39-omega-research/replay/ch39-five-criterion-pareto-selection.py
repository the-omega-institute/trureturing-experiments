"""Joint Pareto, linear-weight and worst-rank selection for the five zero-selection criteria.

On the necessary window 10815..12092 (1278 zeros; 10815..10854 are certified
non-matchable in the 1000 ppm gauge box of the floating-scale certificate) each of
the five published selection criteria is evaluated as a function of the zero label:

    A  four-tap bare spectral weight W(n), maximised (floating-scale winner 10875,
       originally compared only on the gauge-matchable labels),
    B  matching-volume score V(n) = W(n) C(n), maximised (11737; adjustable scale,
       log-uniform top Yukawa, 256 top bins as in ch39-volume-contenders),
    C  the same W up to the constant 13/sqrt(113), originally compared only on
       the six mass candidates 11485..11490 (11488),
    D  least oriented port gain eta(n), minimised (11431); variant: largest
       quasiparticle gap g(n) (also 11431); the finite-temperature parity free
       energy difference is not used as a coordinate,
    E  many-body vacuum energy E(n) of the 38-mode shell, minimised (11618).

Extending A and C beyond their original comparison sets is a model choice and is
stated as such.  All scores are oriented "larger is better" as
s = (W, V, -eta, -E, g) and enclosed by outward-rounded binary64 intervals of the
Arb enclosures; every comparison below is a strict comparison of endpoints, and
every linear certificate is checked in exact rational arithmetic.
"""
import itertools
import math
import time
from fractions import Fraction
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import acb, arb, ctx
from scipy.optimize import linprog

root = Path(__file__).parent
FIRST, COUNT = 10815, 1278
LABELS = list(range(FIRST, FIRST + COUNT))
COLS = ("W", "V", "eta", "E", "g")
MAIN, FIVE = (0, 1, 2, 3), (0, 1, 2, 3, 4)
OMEGA_A = list(range(10855, FIRST + COUNT))      # drop the labels <= 10854 excluded for A
SIX = list(range(11485, 11491))
PIECES = 256


def outward(x):
    lo = math.nextafter(float(x.lower()), -math.inf)
    hi = math.nextafter(float(x.upper()), math.inf)
    assert lo <= hi
    return lo, hi


def collect():
    data = {n: [None] * 5 for n in LABELS}
    vol = run_path(str(root / "ch39-volume-contenders.py"))
    ctx.prec = 112
    gauge_lo, gauge_hi, logs = vol["input_bounds"]()
    prep = vol["prepare"](gauge_lo, gauge_hi, logs)
    targets, _ = vol["zero_targets"](0)
    assert sorted(targets) == LABELS
    q = vol["q"]
    y0 = 2 * (vol["gauge"]["yukawa_centers"][2] * q(".99")).log()
    y1 = 2 * (vol["gauge"]["yukawa_centers"][2] * q("1.01")).log()
    edges = [y0 + (y1 - y0) * j / PIECES for j in range(PIECES + 1)]
    bins = [a.union(b) - logs[2].mid() for a, b in zip(edges, edges[1:])]
    linear = [(b / 2).exp() for b in edges]
    masses = [[arb(1) / PIECES] * PIECES,
              [(b - a) / (linear[-1] - linear[0]) for a, b in zip(linear, linear[1:])]]
    protocols = {}
    for n in LABELS:
        target, weight = targets[n]
        data[n][0] = outward(weight)
        intervals = [vol["widths"](target, prep, y) for y in bins]
        for mode in range(2):
            for measure, mass in enumerate(masses):
                lo = weight * sum(p * v[mode][0] for p, v in zip(mass, intervals))
                hi = weight * sum(p * v[mode][1] for p, v in zip(mass, intervals))
                lo, hi = lo.lower(), hi.upper()
                assert lo >= 0 and lo <= hi
                protocols.setdefault((mode, measure), {})[n] = (
                    math.nextafter(float(lo), -math.inf) if lo > 0 else 0.0,
                    math.nextafter(float(hi), math.inf) if hi > 0 else 0.0)
    for n in LABELS:
        data[n][1] = protocols[0, 0][n]          # adjustable scale, log-uniform top Yukawa
    # The six-candidate weight of ch39-prime-window-six-candidates is (13/sqrt(113)) W.
    ctx.prec = 192
    k = arb(4) / 9
    h, log_two = arb(2).sqrt(), arb(2).log()
    hv = 2 * h * k / (1 + 2 * k * k)
    for gamma in (z.imag for z in acb.zeta_zeros(11485, 6)):
        v = (gamma * log_two).cos()
        six = (((3 * h - 4 * v) / (gamma**2 + arb(1) / 4))**2 * (1 + hv * v)).sqrt()
        u = h * v
        w = (3 - 2 * u) * (2 * (113 + 72 * u)).sqrt() / 13 / (gamma**2 + arb(1) / 4)
        assert abs(six / w - 13 / arb(113).sqrt()) < arb(10)**-50
    gain = run_path(str(root / "ch39-shell-ckm-schur-gain.py"))
    ctx.prec = 128
    near = gain["signed_branches"]()
    for n, z in enumerate(acb.zeta_zeros(FIRST, COUNT), FIRST):
        lo, hi = outward(gain["min_gain"](gain["shell_coefficient"](z.imag), near))
        data[n][2] = (-hi, -lo)
    fan = run_path(str(root / "ch39-parity-weighted-vacuum-fan.py"))
    spectra = fan["collect"]()
    reference = arb(spectra[11618][0].mid())
    sectors = {}
    for n, (energy, g, sector) in spectra.items():
        lo, hi = outward(energy - reference)         # E is shifted by a common exact constant
        data[n][3] = (-hi, -lo)
        data[n][4] = outward(g)
        sectors[n] = sector
    return data, protocols, sectors


def arrays(data, labels):
    lo = np.array([[data[n][j][0] for j in range(5)] for n in labels])
    hi = np.array([[data[n][j][1] for j in range(5)] for n in labels])
    assert np.all(lo <= hi)
    return lo, hi


def pareto(lo, hi, cols):
    """Return per-row status: certainly dominated, certainly optimal, or undecided."""
    lo, hi = lo[:, cols], hi[:, cols]
    status = []
    for k in range(len(lo)):
        weak = np.all(lo >= hi[k], axis=1) & np.any(lo > hi[k], axis=1)
        weak[k] = False
        if weak.any():
            status.append("dominated")
            continue
        possible = np.all(hi >= lo[k], axis=1)
        possible[k] = False
        status.append("undecided" if possible.any() else "optimal")
    return status


def ranks(lo, hi, cols):
    """Worst rank over the criteria, as an interval [certain, possible]."""
    out = []
    for k in range(len(lo)):
        a = [1 + int(np.sum(lo[:, j] > hi[k, j])) for j in cols]
        b = [1 + int(np.sum(np.delete(hi[:, j], k) > lo[k, j])) for j in cols]
        out.append((max(a), max(b), a, b))
    return out


def fractions(lo, hi):
    return ([[Fraction(x) for x in row] for row in lo], [[Fraction(x) for x in row] for row in hi])


def selectability(lo, hi, cols, k, exact):
    """Linear unique-winner test for row k: ('win', w), ('lose', lambda) or ('undecided',)."""
    flo, fhi = exact
    cols = list(cols)
    others = [m for m in range(len(lo)) if m != k]
    base = lo[:, cols].min(0)
    span = hi[:, cols].max(0) - base
    nl, nh = (lo[:, cols] - base) / span, (hi[:, cols] - base) / span
    d = len(cols)
    gap = nl[k][None, :] - nh[others]
    res = linprog(np.r_[np.zeros(d), -1.0], A_ub=np.c_[-gap, np.ones(len(others))],
                  b_ub=np.zeros(len(others)), A_eq=[np.r_[np.ones(d), 0]], b_eq=[1],
                  bounds=[(0, None)] * d + [(None, None)], method="highs")
    if res.status == 0 and -res.fun > 1e-10:
        w = [Fraction(max(x, 0.0)).limit_denominator(10**12) / Fraction(s)
             for x, s in zip(res.x[:d], span)]
        if all(sum(wj * (flo[k][j] - fhi[m][j]) for wj, j in zip(w, cols)) > 0 for m in others):
            return ("win", w)
    res = linprog(np.r_[np.zeros(len(others)), -1.0], A_ub=np.c_[-nl[others].T, np.ones(d)],
                  b_ub=-nh[k], A_eq=[np.r_[np.ones(len(others)), 0]], b_eq=[1],
                  bounds=[(0, None)] * len(others) + [(None, None)], method="highs")
    if res.status == 0 and -res.fun > 1e-10:
        lam = {others[i]: Fraction(x).limit_denominator(10**12)
               for i, x in enumerate(res.x[:-1]) if x > 1e-13}
        total = sum(lam.values())
        lam = {m: x / total for m, x in lam.items()}
        # A convex combination of other rows weakly dominates row k in every column,
        # so no non-negative weight on any subset of these columns makes k a unique winner.
        if all(sum(x * flo[m][j] for m, x in lam.items()) >= fhi[k][j] for j in cols):
            return ("lose", lam)
    return ("undecided",)


def threshold(lo, hi, k, i, j):
    """Bracket [t_lo, t_hi] of the largest t with row k the unique maximiser of s_i + t s_j."""
    flo = lambda m, c: Fraction(lo[m, c])
    fhi = lambda m, c: Fraction(hi[m, c])
    t_lo = t_hi = None
    for m in range(len(lo)):
        if m == k:
            continue
        assert lo[k, i] > hi[m, i]
        if hi[m, j] > lo[k, j]:
            v = (flo(k, i) - fhi(m, i)) / (fhi(m, j) - flo(k, j))
            t_lo = v if t_lo is None or v < t_lo else t_lo
        if lo[m, j] > hi[k, j]:
            v = (fhi(k, i) - flo(m, i)) / (flo(m, j) - fhi(k, j))
            if t_hi is None or v < t_hi:
                t_hi, block = v, m
    assert 0 < t_lo <= t_hi
    return t_lo, t_hi, block


def verify():
    start = time.time()
    data, protocols, sectors = collect()
    print("data collected in", round(time.time() - start), "s")
    lo, hi = arrays(data, LABELS)
    row = {n: k for k, n in enumerate(LABELS)}

    def champion(j, labels):
        rows = [row[n] for n in labels]
        best = max(rows, key=lambda r: lo[r, j])
        assert all(lo[best, j] > hi[r, j] for r in rows if r != best)
        return LABELS[best]

    # 1. Champions on the full window and on the original comparison sets.
    assert [champion(j, LABELS) for j in range(5)] == [10854, 11737, 11431, 11618, 11431]
    assert champion(0, OMEGA_A) == 10875 and champion(0, SIX) == 11488
    assert [champion(j, SIX) for j in range(5)] == [11488, 11488, 11490, 11490, 11490]
    second = sorted(range(COUNT), key=lambda r: -lo[r, 0])[:2]
    assert [LABELS[r] for r in second] == [10854, 10875]
    assert all(lo[row[10875], 0] > hi[r, 0] for r in range(COUNT) if LABELS[r] not in (10854, 10875))
    zero_v = [n for n in LABELS if hi[row[n], 1] == 0]
    assert zero_v == list(range(10815, 10868)) + list(range(12071, 12093))
    assert lo[row[10875], 1] > 0
    print("champions W,V,-eta,-E,g:", [champion(j, LABELS) for j in range(5)],
          "| W on 10855..12092:", champion(0, OMEGA_A), "| W on six:", champion(0, SIX))

    # 2. Pareto classification.
    status = dict(zip(LABELS, pareto(lo, hi, MAIN)))
    counts = {s: sum(v == s for v in status.values()) for s in ("dominated", "optimal", "undecided")}
    undecided = [n for n, s in status.items() if s == "undecided"]
    assert counts == {"dominated": 1168, "optimal": 102, "undecided": 8}
    assert undecided == [10871, 10903, 10972, 11391, 11596, 11722, 11733, 11744]
    for c in (10854, 10875, 11737, 11431, 11618):
        assert status[c] == "optimal"
    assert status[11488] == "dominated"
    strict = [LABELS[m] for m in range(COUNT) if np.all(lo[m] > hi[row[11488]])]
    assert strict == [11575, 11618, 11704]
    # No other label can possibly be strictly better in all five scores.
    assert [m for m in range(COUNT) if m != row[11488] and np.all(hi[m] > lo[row[11488]])] == \
        [row[n] for n in strict]
    # 11618 shares the negative Pfaffian sector of 11488, so the exclusion holds for
    # sector-dependent rules that are non-decreasing in every score (g rewarded).
    assert sectors[11488] == sectors[11618] == -1
    # With a large gap penalised instead, no label of the same sector dominates 11488.
    k8 = row[11488]
    flipped = [m for m in range(COUNT) if m != k8 and sectors[LABELS[m]] == sectors[11488]
               and np.all(hi[m, :4] >= lo[k8, :4]) and -lo[m, 4] >= -hi[k8, 4]]
    assert flipped == []
    # The eight undecided labels are undecided only through the interval of V:
    # every possible dominator is certainly better in each of the other three scores.
    for n in undecided:
        k = row[n]
        rivals = [m for m in range(COUNT) if m != k and np.all(hi[m, :4] >= lo[k, :4])]
        assert rivals and all(np.all(lo[m, [0, 2, 3]] > hi[k, [0, 2, 3]]) for m in rivals)
    status5 = dict(zip(LABELS, pareto(lo, hi, FIVE)))
    c5 = {s: sum(v == s for v in status5.values()) for s in ("dominated", "optimal", "undecided")}
    assert c5 == {"dominated": 1165, "optimal": 104, "undecided": 9}
    loa, hia = arrays(data, OMEGA_A)
    statusa = dict(zip(OMEGA_A, pareto(loa, hia, MAIN)))
    ca = {s: sum(v == s for v in statusa.values()) for s in ("dominated", "optimal", "undecided")}
    assert ca == {"dominated": 1133, "optimal": 98, "undecided": 7}
    assert statusa[11488] == "dominated" and all(statusa[c] == "optimal" for c in (10875, 11737, 11431, 11618))
    los, his = arrays(data, SIX)
    statuss = dict(zip(SIX, pareto(los, his, MAIN)))
    assert statuss == dict(zip(SIX, pareto(los, his, FIVE)))
    assert [n for n in SIX if statuss[n] == "optimal"] == [11488, 11489, 11490]
    assert [n for n in SIX if statuss[n] == "dominated"] == [11485, 11486, 11487]
    print("Pareto (W,V,eta,E) on 1278:", counts, "undecided", undecided)
    print("  with g:", c5, "| on 10855..12092:", ca, "| six optimal: 11488 11489 11490")
    print("  strict dominators of 11488 in all five scores:", strict)

    # Every protocol of the volume score keeps the champion statuses.
    for key, table in protocols.items():
        alt = {n: list(v) for n, v in data.items()}
        for n in LABELS:
            alt[n][1] = table[n]
        alo, ahi = arrays(alt, LABELS)
        st = dict(zip(LABELS, pareto(alo, ahi, MAIN)))
        assert st[11488] == "dominated" and all(st[c] == "optimal" for c in (10875, 11737, 11431, 11618))
        assert all(np.all(alo[row[m]] > ahi[row[11488]]) for m in (11575, 11618, 11704))

    # 3. Linear unique winners.
    exact = fractions(lo, hi)
    result = {}
    for cols, name, expect in ((MAIN, "W,V,eta,E", (29, 23)), (FIVE, "W,V,eta,E,g", (34, 23))):
        st = status if cols == MAIN else status5
        wins, open_ = [], []
        for n in LABELS:
            if st[n] == "dominated":
                continue
            r = selectability(lo, hi, cols, row[n], exact)
            result[cols, n] = r
            (wins if r[0] == "win" else open_ if r[0] == "undecided" else []).append(n)
        assert (len(wins), len(open_)) == expect
        assert sum(st[n] != "dominated" for n in LABELS) - len(wins) - len(open_) == (58 if cols == MAIN else 1278 - 1165 - 34 - 23)
        print("linear unique winners for", name, ":", len(wins), "certified", wins)
        print("  undecided", open_)
        if cols == MAIN:
            main_wins, main_open = wins, open_
            assert wins == [10821, 10839, 10854, 10908, 11032, 11033, 11075, 11176, 11219, 11236,
                            11304, 11305, 11344, 11431, 11491, 11510, 11571, 11575, 11576, 11592,
                            11618, 11619, 11703, 11704, 11726, 11737, 11748, 11891, 11977]
            assert open_ == [10918, 10946, 10968, 11015, 11069, 11172, 11193, 11273, 11290, 11301,
                             11322, 11355, 11381, 11387, 11391, 11398, 11409, 11499, 11521, 11527,
                             11596, 11672, 11744]
    assert open_ == main_open
    # The same classification on Omega_A, and with W replaced by log W.
    exact_a = fractions(loa, hia)
    rowa = {n: k for k, n in enumerate(OMEGA_A)}
    wins_a, open_a = [], []
    for n in OMEGA_A:
        if statusa[n] == "dominated":
            continue
        r = selectability(loa, hia, MAIN, rowa[n], exact_a)
        (wins_a if r[0] == "win" else open_a if r[0] == "undecided" else []).append(n)
    assert wins_a == sorted(set(main_wins) - {10821, 10839, 10854} | {10875}) and len(wins_a) == 27
    assert open_a == main_open
    ctx.prec = 128
    log_lo, log_hi = lo.copy(), hi.copy()
    for k in range(COUNT):
        log_lo[k, 0] = math.nextafter(float(arb(lo[k, 0]).log().lower()), -math.inf)
        log_hi[k, 0] = math.nextafter(float(arb(hi[k, 0]).log().upper()), math.inf)
    assert pareto(log_lo, log_hi, MAIN) == [status[n] for n in LABELS]
    exact_log = fractions(log_lo, log_hi)
    wins_l, open_l = [], []
    for n in LABELS:
        if status[n] == "dominated":
            continue
        r = selectability(log_lo, log_hi, MAIN, row[n], exact_log)
        (wins_l if r[0] == "win" else open_l if r[0] == "undecided" else []).append(n)
    assert (wins_l, open_l) == (main_wins, main_open)
    print("  on 10855..12092:", len(wins_a), "certified (10875 included), same 23 undecided;"
          " log W instead of W: identical lists")
    # Exact range of convex weights x on 10908 with (1-x) s(10854) + x s(10908) >= s(10875).
    a_, b_, c_ = row[10854], row[10908], row[10875]
    low, high = Fraction(0), Fraction(1)
    for j in FIVE:
        p0, p1, need = Fraction(lo[a_, j]), Fraction(lo[b_, j]), Fraction(hi[c_, j])
        if p1 > p0:
            low = max(low, (need - p0) / (p1 - p0))
        elif p1 < p0:
            high = min(high, (p0 - need) / (p0 - p1))
        else:
            assert p0 >= need
    assert Fraction("0.1295") < low < Fraction("0.1296") and Fraction("0.2380") < high < Fraction("0.2381")
    print("  feasible weight on 10908:", float(low), float(high))
    lam = result[FIVE, 10875]
    assert lam[0] == "lose" and set(lam[1]) == {row[10854], row[10908]}
    share = lam[1][row[10908]]
    print("  10875 is weakly dominated by", {LABELS[m]: float(x) for m, x in lam[1].items()})
    assert Fraction(1, 5) < share < Fraction(21, 100) and low <= share <= high
    assert result[FIVE, 10875] and result[MAIN, 10875][0] == "lose"
    for c in (11737, 11431, 11618):
        assert result[MAIN, c][0] == "win" and result[FIVE, c][0] == "win"

    table = []
    for c, i in ((11737, 1), (11431, 2), (11618, 3)):
        for j in MAIN:
            if j != i:
                t_lo, t_hi, block = threshold(lo, hi, row[c], i, j)
                table.append((c, COLS[i], COLS[j], t_lo, t_hi, LABELS[block]))
                sign = lambda x: ("-" if x in ("eta", "E") else "+") + x
                print(f"  {c}: {sign(COLS[i])} {sign(COLS[j])[0]} t*{COLS[j]} unique for t < {float(t_lo):.10g},"
                      f" certainly beaten by {LABELS[block]} for t > {float(t_hi):.10g}")
    expected = {
        (11737, "W"): ("0.00294985", "0.0413077", 11726),
        (11737, "eta"): ("4.99548e-14", "6.99531e-13", 11726),
        (11737, "E"): ("1.47578e-10", "6.63823e-10", 11748),
        (11431, "W"): ("2164989.4", "2164989.5", 11032),
        (11431, "V"): ("207419309", "214549154", 11704),
        (11431, "E"): ("6.48701944", "6.48701945", 11033),
        (11618, "W"): ("176040.80", "176040.82", 11219),
        (11618, "V"): ("56173862", "63126164", 11748),
        (11618, "eta"): ("0.00114697759", "0.0011469776", 11619),
    }
    for c, _, j, t_lo, t_hi, block in table:
        a, b, m = expected[c, j]
        assert Fraction(a) < t_lo <= t_hi < Fraction(b) and block == m, (c, j, float(t_lo), float(t_hi))

    # Cross-check with the diagonal of the parity-weighted fan: 11618 yields to 11619.
    t_lo, t_hi, block = threshold(lo, hi, row[11618], 3, 4)
    assert block == 11619 - FIRST and abs(t_lo - Fraction("0.0673258848")) < Fraction(1, 10**9)
    assert abs(t_hi - Fraction("0.0673258848")) < Fraction(1, 10**9)

    # 4. Worst rank.
    def minimax(lo_, hi_, labels, cols):
        r = ranks(lo_, hi_, cols)
        best = min(range(len(r)), key=lambda k: r[k][1])
        assert all(r[best][1] < r[k][0] for k in range(len(r)) if k != best)
        runner = min((k for k in range(len(r)) if k != best), key=lambda k: r[k][0])
        return labels[best], r[best], labels[runner], r[runner]
    every_lo = [x[0] for x in ranks(lo, hi, MAIN)]
    w, rb, ru, rr = minimax(lo, hi, LABELS, MAIN)
    assert (w, rb[0], rb[1], ru, rr[0], rr[1]) == (11618, 175, 175, 11575, 204, 204)
    third = sorted(range(COUNT), key=lambda k: every_lo[k])[2]
    assert LABELS[third] == 11704 and ranks(lo, hi, MAIN)[third][:2] == (235, 235)
    assert rb[2] == [175, 36, 171, 1] and rb[3] == [175, 39, 171, 1]
    every = ranks(lo, hi, MAIN)
    table_ranks = {10854: ([1, 1195, 424, 1136], [1, 1204, 424, 1136]),
                   10875: ([2, 1181, 575, 1172], [2, 1200, 575, 1172]),
                   11488: ([252, 137, 239, 194], [252, 145, 239, 194]),
                   11737: ([209, 1, 1078, 426], [209, 1, 1078, 426]),
                   11431: ([442, 246, 1, 173], [442, 255, 1, 173])}
    for c, (a_, b_) in table_ranks.items():
        assert every[row[c]][2] == a_ and every[row[c]][3] == b_, c
    assert minimax(lo, hi, LABELS, FIVE)[:2] == (11618, (175, 175, [175, 36, 171, 1, 158], [175, 39, 171, 1, 158]))
    wa = minimax(loa, hia, OMEGA_A, MAIN)
    assert (wa[0], wa[1][:2], wa[2], wa[3][:2]) == (11618, (163, 163), 11575, (191, 191))
    ra = ranks(loa, hia, MAIN)
    thirda = sorted(range(len(OMEGA_A)), key=lambda k: ra[k][0])[2]
    assert OMEGA_A[thirda] == 11704 and ra[thirda][:2] == (220, 220)
    ws = minimax(los, his, SIX, MAIN)
    assert (ws[0], ws[1][:2]) == (11489, (2, 2)) and ws[1][2] == [2, 2, 2, 2]
    rs = ranks(los, his, MAIN)
    assert [x[0] for x in rs] == [x[1] for x in rs] == [6, 6, 4, 3, 2, 5]
    drop = {}
    for cols in itertools.combinations(range(5), 3):
        drop[cols] = minimax(lo, hi, LABELS, cols)[:2]
    assert drop[(0, 2, 3)][0] == 11305 and drop[(0, 2, 3)][1][:2] == (89, 89)
    assert drop[(0, 1, 3)][0] == 11592 and drop[(0, 1, 2)][0] == 11434
    assert drop[(1, 2, 3)][0] == 11575
    for p in protocols:
        alt = {n: list(v) for n, v in data.items()}
        for n in LABELS:
            alt[n][1] = protocols[p][n]
        alo, ahi = arrays(alt, LABELS)
        best_alt = minimax(alo, ahi, LABELS, MAIN)
        assert best_alt[0] == 11618 and best_alt[1][:2] == (175, 175)
    print("worst-rank winner (W,V,eta,E) on 1278: 11618, worst rank 175; runner-up 11575 at 204")
    print("  on 10855..12092: 11618 at 163, runner-up 11575 at 191; on six: 11489 at 2 (11488 at 3)")
    print("  without V: 11305 (89); without eta: 11592; without E: 11434; without W: 11575")
    print("total", round(time.time() - start), "s")


if __name__ == "__main__":
    verify()
