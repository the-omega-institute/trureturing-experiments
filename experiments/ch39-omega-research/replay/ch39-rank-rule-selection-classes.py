"""Which labels each rank-based selection rule can pick: worst rank, Borda and lexicographic.

The five selection readings s = (W, V, -eta, -E, g) and their certified intervals are
those of ch39-five-criterion-pareto-selection, with V intersected as in
ch39-joint-selection-volume-top-split (ladder levels, bisected sub-boxes and the eight
top-Yukawa sub-boxes); the intervals are read from the data files of those two scripts.

Ranks are competition ranks inside the population X (Omega = 10815..12092 or
Omega_A = 10855..12092): r_j(n) = 1 + #{m in X : s_j(m) > s_j(n)}.  An interval
reading gives the rank interval [1 + #{m : lo_j(m) > hi_j(n)}, 1 + #{m != n : hi_j(m) > lo_j(n)}].

Rule classes (w >= 0, w != 0; n is selectable when it is the unique minimiser for some w):

  weighted worst rank   max_j w_j r_j(n);
  weighted Borda        sum_j w_j r_j(n);
  lexicographic         an order of the criteria, ties on one criterion broken by the next.

Certificates:
  * worst rank: the set is the set of labels that no other label weakly dominates (a
    written proof, see the note); with pairwise disjoint W intervals it is the Pareto
    set, whose certified classification is recomputed here;
  * Borda: the linear test of ch39-five-criterion-pareto-selection applied to the rank
    intervals (-r_hi, -r_lo); positive and Farkas certificates are checked in exact
    rational arithmetic; labels that are certainly dominated cannot win any monotone
    rule;
  * the undecided Borda labels: case split on the order of overlapping V pairs; each
    leaf is again an exact certificate;
  * lexicographic: the unique certified champion of the first criterion.

Default run: about 10 s of single-core CPU time, no Arb recomputation of V.
  --width-probe   numerical probe (not a certificate) of the V enclosure on top sub-box 3:
                  variational steps 4096 -> 65536 and two Picard refinements, about 80 s
                  of single-core CPU time.
"""
import itertools
import json
import sys
import time
from fractions import Fraction
from pathlib import Path
from runpy import run_path

import numpy as np
from scipy.optimize import linprog

root = Path(__file__).parent
TS = run_path(str(root / "ch39-joint-selection-volume-top-split.py"))
R2 = TS["R2"]
LAD = TS["LAD"]
OMEGA = list(range(10815, 12093))
OMEGA_A = list(range(10855, 12093))
FOUR, FIVE = [0, 1, 2, 3], [0, 1, 2, 3, 4]
NAMES = ("W", "V", "-eta", "-E", "g")
LINEAR_WINS = sorted(LAD["WINS_MAIN"] + [11290, 11301, 11387, 11672])
LINEAR_UND = [10968, 11069, 11193, 11391, 11499, 11521, 11527, 11596]
BORDA_WINS = [10821, 10839, 10854, 10908, 10946, 10968, 11015, 11032, 11069, 11075, 11107,
              11172, 11176, 11193, 11219, 11236, 11290, 11301, 11303, 11304, 11305, 11344,
              11355, 11381, 11387, 11431, 11491, 11506, 11510, 11571, 11575, 11576, 11592,
              11618, 11619, 11672, 11703, 11704, 11726, 11737, 11747, 11748, 11891, 11977]
BORDA_UND = [11409, 11521, 11527, 11596, 11603, 11640]


def intervals():
    """Certified intervals: published table, ladder and sub-box unions, top-split scores."""
    raw = json.loads(TS["DATA"].read_text())
    new = {int(n): tuple(float.fromhex(x) for x in v) for n, v in raw["scores"].items()}
    table, alt = TS["baseline"]()
    for n, (a, b) in new.items():
        c, d = alt[n][1]
        assert max(a, c) <= min(b, d), n
        alt[n][1] = (max(a, c), min(b, d))
    return alt


def rank_counts(lo, hi, cols):
    """Certain and possible numbers of labels strictly better, per criterion."""
    n = len(lo)
    cert = np.zeros((n, len(cols)), dtype=np.int64)
    poss = np.zeros((n, len(cols)), dtype=np.int64)
    for c, j in enumerate(cols):
        for k in range(n):
            cert[k, c] = np.sum(lo[:, j] > hi[k, j])
            poss[k, c] = np.sum(np.delete(hi[:, j], k) > lo[k, j])
    return cert, poss


def borda_test(rl, rh, k):
    """('win', w) with sum w r_hi(k) < sum w r_lo(m) for all m != k; ('lose', lam) with a
    convex combination of the r_hi(m) at most r_lo(k) in every criterion; else undecided."""
    n, d = rl.shape
    oth = [m for m in range(n) if m != k]
    A = rl[oth] - rh[k]
    scale = np.maximum(np.abs(A).max(0), 1)
    res = linprog(np.r_[np.zeros(d), -1.0], A_ub=np.c_[-A / scale, np.ones(len(oth))],
                  b_ub=np.zeros(len(oth)), A_eq=[np.r_[np.ones(d), 0]], b_eq=[1],
                  bounds=[(0, None)] * d + [(None, None)], method="highs")
    if res.status == 0 and -res.fun > 1e-12:
        w = [Fraction(max(x, 0.0)).limit_denominator(10**9) / int(s) for x, s in zip(res.x[:d], scale)]
        if all(sum(wj * int(a) for wj, a in zip(w, row)) > 0 for row in A):
            return "win", w
    B = rh[oth] - rl[k]
    res = linprog(np.r_[np.zeros(len(oth)), 1.0], A_ub=np.c_[B.T, -np.ones(d)], b_ub=np.zeros(d),
                  A_eq=[np.r_[np.ones(len(oth)), 0]], b_eq=[1],
                  bounds=[(0, None)] * len(oth) + [(None, None)], method="highs")
    if res.status == 0 and res.fun < 1e-9:
        lam = {oth[i]: Fraction(x).limit_denominator(10**9) for i, x in enumerate(res.x[:-1]) if x > 1e-12}
        total = sum(lam.values())
        lam = {m: x / total for m, x in lam.items()}
        if all(sum(x * int(rh[m, c] - rl[k, c]) for m, x in lam.items()) <= 0 for c in range(d)):
            return "lose", lam
    return "undecided", None


def classify(alt, labels, cols):
    lo, hi = R2["arrays"](alt, labels)
    pareto = R2["pareto"](lo, hi, cols)
    cert, poss = rank_counts(lo, hi, cols)
    out = {}
    for k, n in enumerate(labels):
        # A label that is certainly weakly dominated is never the unique minimiser of a
        # rule that is monotone in every rank: its dominator has no larger rank anywhere.
        out[n] = ("lose", "dominated") if pareto[k] == "dominated" else borda_test(1 + cert, 1 + poss, k)
    return lo, hi, pareto, cert, poss, out


def overlapping_pairs(lo, hi, j):
    """Pairs whose order in criterion j is undecided by the intervals."""
    out = []
    for a in range(len(lo)):
        for b in range(a + 1, len(lo)):
            if hi[a, j] > lo[b, j] and hi[b, j] > lo[a, j] and not (lo[a, j] > hi[b, j] or lo[b, j] > hi[a, j]):
                out.append((a, b))
    return out


def split(k, cert, poss, c, pairs, order):
    """Borda test of row k after fixing the V order of the given pairs (case split).

    order[(a, b)] is '>' (V(a) > V(b)), '=' or '<'.  Fixing a > b makes a certainly
    better than b and removes b from the possible betters of a; '=' removes each from
    the possible betters of the other.  Every true configuration lies in one case."""
    C, P = cert.copy(), poss.copy()
    for (a, b), o in order.items():
        assert (a, b) in pairs
        if o == ">":
            C[b, c] += 1
            P[a, c] -= 1
        elif o == "<":
            C[a, c] += 1
            P[b, c] -= 1
        else:
            P[a, c] -= 1
            P[b, c] -= 1
    return borda_test(1 + C, 1 + P, k)[0]


def width_probe():
    """Numerical probe only: V enclosure on top sub-box 3 with more variational steps."""
    from flint import arb, arb_mat
    vol, gauge = TS["vol"], TS["vol"]["gauge"]
    M = TS["Model"]()
    src = (root / "ch39-volume-contenders.py").read_text()
    code = src[src.index("def variations"):src.index("def prepare")]
    assert code.count("anchor+radius, 4096,") == 1
    k = 3
    e0, e1 = TS["sub_edges"](M, k)
    logs = list(M.logs)
    logs[2] = (e0 - TS["EPS"]).union(e1 + TS["EPS"])
    gl, gh = M.gl, M.gh
    low = [a.union(b) for a, b in zip(gl, gh)]
    xr = [(a + b) / 2 for a, b in zip(gl, gh)]
    yr = [y.mid() for y in logs]
    end, _, _ = gauge["evolve"](xr, yr, vol["anchor"], TS["CENTER_STEPS"], diagonal_quarks=True)
    raw = json.loads(TS["DATA"].read_text())
    ref0, P0 = TS["des_bounds"](raw["preps"][k])
    ref = logs[2].mid()

    def width(P, r, n):
        target, _ = M.target(n)
        a = b = arb(0)
        edges = [e0 + (e1 - e0) * j / 128 for j in range(129)]
        for x0, x1 in zip(edges, edges[1:]):
            l, u = LAD["parts"](target, P, x0.union(x1) - r, vol["radius"], 16)
            a, b = a + l, b + u
        return float((b - a) / b)

    out = {}
    for steps, refinements in ((4096, 1), (65536, 2)):
        ns = dict(vol)
        exec(code.replace("anchor+radius, 4096,", f"anchor+radius, {steps},"), ns)
        jac, rates = ns["variations"](low, logs, refinements=refinements)
        m = arb_mat([[jac[i, j] for j in range(3)] for i in range(3)])
        pre = arb_mat([[m[i, j].mid() for j in range(3)] for i in range(3)]).inv()
        det = m.det().intersection((pre * m).det() / pre.det())
        col = lambda v: arb_mat([[x] for x in v])
        normalized, slope, center = pre * col(vol["ratios"]), pre * col(rates), pre * col(end)
        rem = pre * jac
        for i in range(3):
            rem[i, i] -= 1
        inputs = [v - x for v, x in zip(low, xr)] + [y - r for y, r in zip(logs, yr)]
        res = [sum(rem[i, j] * inputs[j] for j in range(12) if j != 5) for i in range(3)]
        faces = [(center[i, 0] + res[i] + (gh[i] - xr[i]), center[i, 0] + res[i] + (gl[i] - xr[i]))
                 for i in range(3)]
        P = TS["Bounds"](faces, [rem[i, 5] for i in range(3)], [normalized[i, 0] for i in range(3)],
                         [slope[i, 0] for i in range(3)], rates, det)
        if steps == 4096:
            # Same computation as the stored preparation: same widths up to the decoding pad.
            assert abs(width(P, ref, 11521) - width(P0, ref0, 11521)) < 1e-6
        out[steps] = {"residual": max(float(x.rad()) for x in res),
                      "center": max(float(center[i, 0].rad()) for i in range(3)),
                      "det": float(det.rad()),
                      "w": {n: width(P, ref, n) for n in (11521, 11409)}}
        print("variational steps", steps, out[steps], flush=True)
    a, b = out[4096], out[65536]
    assert 2.8e-7 < a["residual"] < 2.9e-7 and 2.2e-7 < b["residual"] < 2.3e-7
    assert 1.0e-7 < a["center"] < 1.01e-7 and 1.0e-7 < b["center"] < 1.01e-7
    assert 1.3e-4 < a["det"] < 1.4e-4 and 9.7e-5 < b["det"] < 9.9e-5
    assert 0.76 < b["residual"] / a["residual"] < 0.77
    assert 9.4e-4 < a["w"][11521] < 9.42e-4
    for n in (11521, 11409):
        assert 0.81 < b["w"][n] / a["w"][n] < 0.83, (n, b["w"][n] / a["w"][n])
    print("probe: residual radius stays", b["residual"], "; relative width ratio",
          {n: round(b["w"][n] / a["w"][n], 3) for n in b["w"]})


def verify():
    start = time.time()
    alt = intervals()
    summary = {}

    # 1. Pareto sets and pairwise disjoint W intervals (no two identical score vectors).
    for name, labels, cols, expect in (("Omega four", OMEGA, FOUR, (1172, 106)),
                                       ("Omega_A four", OMEGA_A, FOUR, (1137, 101)),
                                       ("Omega five", OMEGA, FIVE, (1170, 108))):
        lo, hi = R2["arrays"](alt, labels)
        st = R2["pareto"](lo, hi, cols)
        assert (st.count("dominated"), st.count("optimal"), st.count("undecided")) == expect + (0,), name
        o = np.argsort(lo[:, 0])
        assert np.all(lo[o[1:], 0] > hi[o[:-1], 0])
        summary[name] = {"pareto": expect[1]}
    print("worst-rank-selectable = Pareto set:", {k: v["pareto"] for k, v in summary.items()})

    # 2. Lexicographic: the unique certified champion of the leading criterion.
    lex = {}
    for name, labels, cols in (("Omega four", OMEGA, FOUR), ("Omega_A four", OMEGA_A, FOUR),
                               ("Omega five", OMEGA, FIVE)):
        lo, hi = R2["arrays"](alt, labels)
        champ = {}
        for j in cols:
            b = int(np.argmax(lo[:, j]))
            assert all(lo[b, j] > hi[m, j] for m in range(len(labels)) if m != b)
            champ[NAMES[j]] = labels[b]
        lex[name] = champ
        wins = sum(champ[NAMES[p[0]]] == 11618 for p in itertools.permutations(cols))
        lex[name, "11618 orders"] = (wins, len(list(itertools.permutations(cols))))
    assert lex["Omega four"] == {"W": 10854, "V": 11737, "-eta": 11431, "-E": 11618}
    assert lex["Omega_A four"] == {"W": 10875, "V": 11737, "-eta": 11431, "-E": 11618}
    assert lex["Omega five"] == {"W": 10854, "V": 11737, "-eta": 11431, "-E": 11618, "g": 11431}
    assert lex["Omega four", "11618 orders"] == (6, 24) and lex["Omega five", "11618 orders"] == (24, 120)
    assert set(lex["Omega four"].values()) <= set(LINEAR_WINS)
    assert alt[10854][1] == (0.0, 0.0) and min(OMEGA_A) == 10855
    print("lexicographic champions:", {k: v for k, v in lex.items()})

    # 3. Weighted Borda.
    borda = {}
    for name, labels, cols in (("Omega four", OMEGA, FOUR), ("Omega_A four", OMEGA_A, FOUR),
                               ("Omega five", OMEGA, FIVE)):
        lo, hi, st, cert, poss, out = classify(alt, labels, cols)
        borda[name] = (lo, hi, cert, poss, out)
        w = [n for n, r in out.items() if r[0] == "win"]
        u = [n for n, r in out.items() if r[0] == "undecided"]
        assert all(st[labels.index(n)] == "optimal" for n in w + u)
        print(name, "Borda: win", len(w), "undecided", u, "lose", len(labels) - len(w) - len(u),
              round(time.time() - start), "s", flush=True)
        summary[name]["borda"] = (w, u)
    w4, u4 = summary["Omega four"]["borda"]
    assert w4 == BORDA_WINS and u4 == BORDA_UND
    wa, ua = summary["Omega_A four"]["borda"]
    assert len(wa) == 43 and ua == [11409, 11527, 11596, 11603, 11640]
    assert set(wa) == set(BORDA_WINS) - {10821, 10839, 10854} | {10875, 11521}
    w5, u5 = summary["Omega five"]["borda"]
    assert len(w5) == 49 and u5 == BORDA_UND
    assert set(w5) == set(BORDA_WINS) | {11033, 11363, 11635, 11763, 12036}

    # Linear and Borda classes on Omega four are incomparable.
    out4 = borda["Omega four"][4]
    assert sorted(set(LINEAR_WINS) - set(w4)) == [11033]
    only_borda = sorted(set(w4) - set(LINEAR_WINS))
    assert only_borda == [10946, 10968, 11015, 11069, 11107, 11172, 11193, 11303, 11355, 11381,
                          11506, 11747]
    lin_lose_borda_win = [n for n in only_borda if n not in LINEAR_UND]
    assert lin_lose_borda_win == [10946, 11015, 11107, 11172, 11303, 11355, 11381, 11506, 11747]
    status8 = {n: out4[n][0] for n in LINEAR_UND}
    assert status8 == {10968: "win", 11069: "win", 11193: "win", 11391: "lose", 11499: "lose",
                       11521: "undecided", 11527: "undecided", 11596: "undecided"}
    assert out4[11033][0] == "lose" and len(set(w4) & set(LINEAR_WINS)) == 32
    print("linear-only", [11033], "Borda-only", only_borda, "Borda status of the 8 linear-undecided", status8)

    # 4. Equal weights.
    lo, hi, cert, poss, _ = borda["Omega four"]
    s_lo, s_hi = (1 + cert).sum(1), (1 + poss).sum(1)
    k = OMEGA.index(11704)
    assert s_lo[k] == s_hi[k] == 356 and all(s_hi[k] < s_lo[m] for m in range(len(OMEGA)) if m != k)
    second = sorted(range(len(OMEGA)), key=lambda m: s_lo[m])[1]
    assert OMEGA[second] == 11575 and (s_lo[second], s_hi[second]) == (361, 362)
    assert all(s_hi[second] < s_lo[m] for m in range(len(OMEGA)) if m not in (k, second))
    k18 = OMEGA.index(11618)
    assert (s_lo[k18], s_hi[k18]) == (384, 385)
    assert s_hi[second] < s_lo[k18]
    assert all(s_hi[k18] < s_lo[m] for m in range(len(OMEGA)) if m not in (k, second, k18))
    assert list(1 + cert[k]) == list(1 + poss[k]) == [235, 11, 8, 102]
    print("equal-weight Borda: 11704 unique with rank sum 356, then 11575 in [361, 362]; 11618 in [384, 385]")
    for name, labels, cols, sums in (("Omega_A four", OMEGA_A, FOUR, (338, 11575, (345, 346))),
                                     ("Omega five", OMEGA, FIVE, (374, 11575, (433, 434)))):
        lo_, hi_, c_, p_, _ = borda[name]
        a_, b_ = (1 + c_).sum(1), (1 + p_).sum(1)
        kk = labels.index(11704)
        assert a_[kk] == b_[kk] == sums[0] and all(b_[kk] < a_[m] for m in range(len(labels)) if m != kk)
        sec = sorted(range(len(labels)), key=lambda m: a_[m])[1]
        assert labels[sec] == sums[1] and (a_[sec], b_[sec]) == sums[2], (name, labels[sec], a_[sec])
        assert all(b_[sec] < a_[m] for m in range(len(labels)) if m not in (kk, sec))
    assert LAD["minimax"](lo, hi, OMEGA, FOUR) == (11618, (175, 175))

    # 5. The six undecided Borda labels: case split on overlapping V pairs.
    c = 1
    pairs = overlapping_pairs(lo, hi, 1)
    assert len(pairs) == 757
    row = {n: i for i, n in enumerate(OMEGA)}
    P = lambda a, b: (row[a], row[b])
    cases = {}

    def tree(n, plan):
        """plan: list of (pairs fixed so far, status) leaves; returns them checked."""
        leaves = []
        for fixed, expect in plan:
            order = {P(a, b): o for a, b, o in fixed}
            got = split(row[n], cert, poss, c, pairs, order)
            assert got == expect, (n, fixed, got)
            leaves.append((fixed, got))
        return leaves

    def full(first, rest):
        """All orders of `first`, each combined with the list `rest[o]` of sub-cases."""
        out = []
        for o in (">", "=", "<"):
            for sub, status in rest[o]:
                out.append(([(first[0], first[1], o)] + sub, status))
        return out

    W, L, U = "win", "lose", "undecided"
    # 11409 wins iff V(11409) >= V(11533) or V(11409) >= V(11548).
    cases[11409] = tree(11409, full((11409, 11533), {
        ">": [([], W)], "=": [([], W)],
        "<": [([(11409, 11548, ">")], W), ([(11409, 11548, "=")], W), ([(11409, 11548, "<")], L)]}))
    # 11521 wins iff V(11521) >= V(11560) or V(11521) >= V(11781).
    cases[11521] = tree(11521, full((11521, 11560), {
        ">": [([], W)], "=": [([], W)],
        "<": [([(11521, 11781, ">")], W), ([(11521, 11781, "=")], W), ([(11521, 11781, "<")], L)]}))
    # 11527 wins iff V(11521) < V(11560) or V(11521) < V(11781).
    sub = [([(11521, 11781, ">")], L), ([(11521, 11781, "=")], L), ([(11521, 11781, "<")], W)]
    cases[11527] = tree(11527, full((11521, 11560), {">": sub, "=": sub, "<": [([], W)]}))
    # 11603 wins iff V(11603) >= V(11624) or V(11603) >= V(11753) or V(11592) < V(11607).
    cases[11603] = tree(11603, full((11603, 11624), {
        ">": [([], W)], "=": [([], W)],
        "<": [([(11603, 11753, ">")], W), ([(11603, 11753, "=")], W),
              ([(11603, 11753, "<"), (11592, 11607, ">")], L),
              ([(11603, 11753, "<"), (11592, 11607, "=")], L),
              ([(11603, 11753, "<"), (11592, 11607, "<")], W)]}))
    # 11640 wins if V(11640) >= V(11678); otherwise undecided by these splits.
    cases[11640] = tree(11640, [([(11640, 11678, ">")], W), ([(11640, 11678, "=")], W),
                                ([(11640, 11678, "<")], U)])
    # 11596 loses if V(11596) < V(11759); otherwise undecided by this split.
    cases[11596] = tree(11596, [([(11596, 11759, ">")], U), ([(11596, 11759, "=")], U),
                                ([(11596, 11759, "<")], L)])
    for n in BORDA_UND:
        involved = [p for p in pairs if row[n] in p]
        print(n, "overlapping V pairs", [(OMEGA[a], OMEGA[b]) for a, b in involved],
              "leaves", [(f, s) for f, s in cases[n]])
    assert {n: sorted(OMEGA[b] if OMEGA[a] == n else OMEGA[a] for a, b in pairs if row[n] in (a, b))
            for n in BORDA_UND} == {11409: [11533, 11548], 11521: [11560, 11781], 11527: [],
                                    11596: [11759], 11603: [11624, 11753], 11640: [11678]}
    rel = {}
    for a, b in ((11409, 11533), (11409, 11548), (11521, 11560), (11521, 11781), (11603, 11624),
                 (11603, 11753), (11592, 11607), (11640, 11678), (11596, 11759)):
        ma, mb = sum(alt[a][1]) / 2, sum(alt[b][1]) / 2
        rel[a, b] = (ma - mb) / ma
    assert 1.8e-5 < abs(rel[11521, 11560]) < 1.9e-5 and all(2.5e-4 < abs(v) < 9e-4 for p, v in rel.items()
                                                 if p != (11521, 11560))
    print("relative midpoint differences of the deciding pairs (numerical, not a certificate):",
          {p: f"{v:+.2e}" for p, v in rel.items()})

    # 6. 11618 is not a winner for every weight in any class: with w = e_W both rank rules
    #    reduce to r_W and the linear rule to W, whose unique certified champion is 10854;
    #    with w = e_E all three select 11618, the unique champion of -E.
    for label, c in ((10854, 0), (11618, 3)):
        kk = OMEGA.index(label)
        assert cert[kk, c] == poss[kk, c] == 0
        assert all(cert[m, c] >= 1 for m in range(len(OMEGA)) if m != kk)
    assert out4[10854][0] == "win" and out4[11618][0] == "win"

    # 7. Numerical probe only (not a certificate): shrink every V interval about its
    #    midpoint and rerun the four-score linear test on the 8 linear-undecided labels.
    exact = R2["fractions"](lo, hi)
    assert all(R2["selectability"](lo, hi, FOUR, row[n], exact)[0] == "undecided" for n in LINEAR_UND)
    probe = {}
    for factor in (0.1, 0.05):
        l2, h2 = lo.copy(), hi.copy()
        mid, half = (lo[:, 1] + hi[:, 1]) / 2, (hi[:, 1] - lo[:, 1]) / 2 * factor
        l2[:, 1], h2[:, 1] = mid - half, mid + half
        ex2 = R2["fractions"](l2, h2)
        probe[factor] = [n for n in LINEAR_UND if R2["selectability"](l2, h2, FOUR, row[n], ex2)[0] == "undecided"]
    assert probe == {0.1: [11069, 11527], 0.05: []}
    print("linear probe, V half-widths x 0.1 leaves", probe[0.1], "undecided; x 0.05 leaves none")
    print("total", round(time.time() - start), "s")


if __name__ == "__main__":
    if "--width-probe" in sys.argv[1:]:
        width_probe()
    else:
        verify()
