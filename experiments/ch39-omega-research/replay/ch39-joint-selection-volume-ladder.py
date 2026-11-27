"""Refinement ladder for the volume score V of the joint five-criterion selection.

The joint comparison in ch39-five-criterion-pareto-selection leaves 8 Pareto labels and
23 linear-weight labels undecided, all through the enclosure of V.  V's enclosure in
ch39-volume-contenders is the time integral of a density over a matching window: the
lower end integrates the density lower bound over the sufficient window, the upper
end the density upper bound over the necessary window.  Both hold uniformly for every
fixed value of the eight non-top Yukawa inputs.

Ladder level l (l = 0..3) uses 256 * 2**l top-Yukawa bins and 2**l time slabs for the
density; the certified interval at level l is the intersection of levels 0..l (each is
a valid enclosure).  Two quantities are laddered:

  full    the original score, uniform over the eight other Yukawa inputs;
  center  the score with the eight other Yukawa inputs fixed at the box center, an
          admissible instance of the original "any fixed value" reading and a stated
          model choice; its value lies inside every enclosure of the full score.

A third refinement bisects each of the eight non-top Yukawa coordinates once (256
sub-boxes, one `prepare` each) for 10903, 10908, 10968 and 10972; the union of the
sub-box enclosures is again an enclosure of V for every fixed value of those inputs.

Default run: loads the data file, recomputes both ladders on the 15 labels that decide
the Pareto questions and four of the 256 sub-boxes, checks bit equality, reruns the
probes, and redoes every classification from the data.  Writing modes (they overwrite
the data file unless `--output PATH` is given):
  --regenerate            table and both ladders (about seven minutes), keeps sub-boxes;
  --subboxes FIRST STOP   sub-boxes FIRST..STOP-1 (about 3.7 s each), no verification.
"""
import itertools
import json
import math
import sys
import time
from fractions import Fraction
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import acb, arb, ctx

root = Path(__file__).parent
DATA = root / "ch39-joint-selection-volume-ladder-data.json"
R2 = run_path(str(root / "ch39-five-criterion-pareto-selection.py"))
vol = run_path(str(root / "ch39-volume-contenders.py"))
LABELS, MAIN, FIVE, OMEGA_A = R2["LABELS"], R2["MAIN"], R2["FIVE"], R2["OMEGA_A"]
LEVELS = 4
PARETO_UND = [10871, 10903, 10972, 11391, 11596, 11722, 11733, 11744]
RIVALS = {10871: [10854], 10903: [10908], 10972: [10968], 11391: [11387], 11596: [11592],
          11722: [11748], 11733: [11726, 11748], 11744: [11722, 11748]}
LINEAR_UND = [10918, 10946, 10968, 11015, 11069, 11172, 11193, 11273, 11290, 11301, 11322,
              11355, 11381, 11387, 11391, 11398, 11409, 11499, 11521, 11527, 11596, 11672, 11744]
WINS_MAIN = [10821, 10839, 10854, 10908, 11032, 11033, 11075, 11176, 11219, 11236, 11304,
             11305, 11344, 11431, 11491, 11510, 11571, 11575, 11576, 11592, 11618, 11619,
             11703, 11704, 11726, 11737, 11748, 11891, 11977]
SLICE = sorted(set(PARETO_UND) | {m for v in RIVALS.values() for m in v})
PAIRS = [(11391, 11387), (11596, 11592), (11722, 11748), (11733, 11726), (11733, 11748),
         (11744, 11748), (11744, 11722), (10972, 10968), (10903, 10908)]
# Relative overlap of the two V enclosures at levels 0 and 3, and its geometric limit.
EXPECT_OVERLAP = {
    "full": {(11391, 11387): (3.15e-3, 2.34e-3, 2.22e-3), (11596, 11592): (3.25e-3, 2.30e-3, 2.17e-3),
             (11733, 11726): (2.28e-3, 1.22e-3, 1.07e-3), (11733, 11748): (4.80e-3, 3.75e-3, 3.60e-3),
             (11744, 11748): (3.77e-3, 2.72e-3, 2.57e-3), (11744, 11722): (2.40e-3, 1.36e-3, 1.21e-3),
             (10972, 10968): (9.09e-3, 7.63e-3, 7.42e-3), (10903, 10908): (1.71e-2, 1.28e-2, 1.22e-2)},
    "center": {(11733, 11748): (3.08e-3, 2.03e-3, 1.88e-3), (11744, 11748): (2.38e-3, 1.34e-3, 1.19e-3)}}


def fbits(x):
    return float.hex(x)


def unbits(s):
    return float.fromhex(s)


def parts(target, data, top, radius, slabs):
    """Sufficient/necessary-window bounds of the density integral (widths of
    ch39-volume-contenders at zero tolerance, adjustable scale), with the density
    bounded on `slabs` equal time slabs instead of on the whole window."""
    ratios = vol["ratios"]
    both = [[], []]
    for i in range(3):
        norm = data.target_ratios[i]
        assert norm > 0
        for side, value in enumerate(data.faces[i]):
            start = (value + data.top_sensitivity[i] * top) / norm
            rate = data.transformed_rates[i] / norm
            assert target * (2 * radius).exp() > start - radius * rate
            assert target * (-2 * radius).exp() < start + radius * rate
            assert rate + 2 * target * (-2 * radius).exp() > 0
            delta = (-radius).union(radius)
            for _ in range(6):
                delta = (target / (start + rate * delta)).log() / 2
            both[side].append(delta)
    early, late = both
    n0, n1 = max(x.lower() for x in early), min(x.upper() for x in late)
    s0, s1 = max(x.upper() for x in early), min(x.lower() for x in late)

    def density(t):
        common = target * (-2 * t).exp()
        speeds = [2 * common + v / r for v, r in zip(data.rates, ratios)]
        assert all(s > 0 for s in speeds)
        return sum(speeds) / data.determinant

    lo, up = arb(0), arb(0)
    if n1 <= n0:
        return lo, up
    step = (n1 - n0) / slabs
    for k in range(slabs):
        up += step * density((n0 + k * step).union(n0 + (k + 1) * step)).upper()
    if s1 > s0:
        step = (s1 - s0) / slabs
        for k in range(slabs):
            lo += step * density((s0 + k * step).union(s0 + (k + 1) * step)).lower()
    return lo, up


def sufficient_margin(target, data, top, radius):
    """Length s1 - s0 of the sufficient window (negative when it is empty)."""
    both = [[], []]
    for i in range(3):
        norm = data.target_ratios[i]
        for side, value in enumerate(data.faces[i]):
            start = (value + data.top_sensitivity[i] * top) / norm
            rate = data.transformed_rates[i] / norm
            delta = (-radius).union(radius)
            for _ in range(6):
                delta = (target / (start + rate * delta)).log() / 2
            both[side].append(delta)
    early, late = both
    s0 = early[0].upper()
    for x in early[1:]:
        s0 = s0.max(x.upper())
    s1 = late[0].lower()
    for x in late[1:]:
        s1 = s1.min(x.lower())
    return s1 - s0


class Model:
    def __init__(self):
        ctx.prec = 112
        q = vol["q"]
        self.q = q
        self.gl, self.gh, self.logs = vol["input_bounds"]()
        fixed = vol["gauge"]["fixed"]
        self.common = ((fixed["geo"]["BOOK_PLANCK"] / 1000)**2 * 3 * arb(5).sqrt()
                       / (4 * arb.pi()) * (-2 * vol["anchor"]).exp())
        self.offset = q(fixed["OFFSET"])
        self.weight = vol["base"]["spectral"]["weight"]
        self.y0 = 2 * (vol["gauge"]["yukawa_centers"][2] * q(".99")).log()
        self.y1 = 2 * (vol["gauge"]["yukawa_centers"][2] * q("1.01")).log()
        self.center_logs = [y if a == 2 else arb(y.mid()) for a, y in enumerate(self.logs)]
        self.preps = {"full": vol["prepare"](self.gl, self.gh, self.logs),
                      "center": vol["prepare"](self.gl, self.gh, self.center_logs)}

    def target(self, n):
        g = acb.zeta_zero(n).imag
        return self.common / (g * g + self.offset), self.weight(g)

    def score(self, n, prep, pieces, slabs, radius=None):
        radius = vol["radius"] if radius is None else radius
        target, w = self.target(n)
        edges = [self.y0 + (self.y1 - self.y0) * j / pieces for j in range(pieces + 1)]
        a, b = arb(0), arb(0)
        for e0, e1 in zip(edges, edges[1:]):
            lo, up = parts(target, prep, e0.union(e1) - self.logs[2].mid(), radius, slabs)
            a, b = a + lo, b + up
        lo, up = (w * a / pieces).lower(), (w * b / pieces).upper()
        return (math.nextafter(float(lo), -math.inf) if lo > 0 else 0.0,
                math.nextafter(float(up), math.inf) if up > 0 else 0.0)

    def ladder(self, labels, mode):
        return {lev: {n: self.score(n, self.preps[mode], 256 * 2**lev, 2**lev) for n in labels}
                for lev in range(LEVELS)}


def regenerate(path):
    """Recompute table and ladders; keep the separately generated sub-box block."""
    data, protocols, sectors = R2["collect"]()
    lo, hi = R2["arrays"](data, LABELS)
    st = R2["pareto"](lo, hi, MAIN)
    st5 = R2["pareto"](lo, hi, FIVE)
    S = sorted({n for n, a, b in zip(LABELS, st, st5) if a != "dominated" or b != "dominated"}
               | set(SLICE))
    model = Model()
    out = {"table": {str(n): [[fbits(x) for x in data[n][j]] for j in range(5)] for n in LABELS},
           "S": S}
    for mode in ("full", "center"):
        lad = model.ladder(S, mode)
        out[mode] = {str(l): {str(n): [fbits(x) for x in v] for n, v in lad[l].items()} for l in lad}
    old = json.loads(path.read_text()) if path.exists() else {}
    if "subboxes" in old:
        out["subboxes"] = old["subboxes"]
    path.write_text(json.dumps(out, sort_keys=True))
    print("wrote", path)


# Eight non-top Yukawa coordinates, each bisected once: 256 sub-boxes.
NONTOP = [0, 1, 3, 4, 5, 6, 7, 8]
SUB_LABELS = [10903, 10908, 10968, 10972]
SUB_PIECES, SUB_SLABS = 1024, 4


def subbox_logs(model, b):
    out = list(model.logs)
    for j, a in enumerate(NONTOP):
        y = model.logs[a]
        lo, hi = arb(y.lower()), arb(y.upper())
        mid = (lo + hi) / 2
        out[a] = lo.union(mid) if (b >> j) & 1 == 0 else mid.union(hi)
    return out


def subbox_scores(model, b):
    prep = vol["prepare"](model.gl, model.gh, subbox_logs(model, b))
    return {n: model.score(n, prep, SUB_PIECES, SUB_SLABS) for n in SUB_LABELS}


def generate_subboxes(path, first, stop):
    model = Model()
    raw = json.loads(path.read_text())
    block = raw.setdefault("subboxes", {})
    for b in range(first, stop):
        block[str(b)] = {str(n): [fbits(x) for x in v] for n, v in subbox_scores(model, b).items()}
        print("sub-box", b, flush=True)
    path.write_text(json.dumps(raw, sort_keys=True))
    print("wrote sub-boxes", first, "to", stop - 1, "into", path)


def load(path=DATA):
    raw = json.loads(path.read_text())
    table = {int(n): [tuple(unbits(x) for x in pair) for pair in v] for n, v in raw["table"].items()}
    ladders = {mode: {int(l): {int(n): tuple(unbits(x) for x in v) for n, v in raw[mode][l].items()}
                      for l in raw[mode]} for mode in ("full", "center")}
    subs = {int(b): {int(n): tuple(unbits(x) for x in v) for n, v in d.items()}
            for b, d in raw.get("subboxes", {}).items()}
    return table, raw["S"], ladders, subs


def refined(table, S, ladder, level):
    """Five-criterion table with V intersected over ladder levels 0..level on S."""
    alt = {n: list(v) for n, v in table.items()}
    for n in S:
        a, b = table[n][1]
        for l in range(level + 1):
            c, d = ladder[l][n]
            a, b = max(a, c), min(b, d)
        assert a <= b, n
        alt[n][1] = (a, b)
    return alt


def classify_linear(lo, hi, cols, labels, status):
    exact = R2["fractions"](lo, hi)
    row = {n: k for k, n in enumerate(labels)}
    res = {}
    for n in labels:
        if status[n] != "dominated":
            res[n] = R2["selectability"](lo, hi, cols, row[n], exact)[0]
    return res


def minimax(lo, hi, labels, cols):
    r = R2["ranks"](lo, hi, cols)
    best = min(range(len(r)), key=lambda k: r[k][1])
    assert all(r[best][1] < r[k][0] for k in range(len(r)) if k != best)
    return labels[best], r[best][:2]


def verify(path=DATA):
    start = time.time()
    table, S, ladders, subs = load(path)
    assert sorted(table) == LABELS and set(SLICE) <= set(S) and len(S) == 113 and len(SLICE) == 15
    lo0, hi0 = R2["arrays"](table, LABELS)
    st0 = dict(zip(LABELS, R2["pareto"](lo0, hi0, MAIN)))
    assert [n for n in LABELS if st0[n] == "undecided"] == PARETO_UND
    assert sum(v == "dominated" for v in st0.values()) == 1168
    row = {n: k for k, n in enumerate(LABELS)}
    for n, rv in RIVALS.items():
        k = row[n]
        poss = [LABELS[m] for m in range(len(LABELS)) if m != k and np.all(hi0[m, :4] >= lo0[k, :4])]
        assert poss == rv
        # Each possible dominator is already certainly better in W, -eta and -E.
        assert all(np.all(lo0[row[m], [0, 2, 3]] > hi0[k, [0, 2, 3]]) for m in rv)

    # 1. Bit-exact recomputation of both ladders on the decisive slice.
    model = Model()
    for mode in ("full", "center"):
        lad = model.ladder(SLICE, mode)
        for l in range(LEVELS):
            for n in SLICE:
                assert lad[l][n] == ladders[mode][l][n], (mode, l, n)
    print("ladders recomputed on", len(SLICE), "labels in", round(time.time() - start), "s")

    # 2. Where the width of V comes from.
    # Level 0 is already slightly tighter than the published 256-bin interval.
    assert 1.9751394e-10 < table[11391][1][0] < 1.9751395e-10
    assert 1.9751466e-10 < ladders["full"][0][11391][0] < 1.9751467e-10
    rel = lambda iv: (iv[1] - iv[0]) / iv[1]
    full = ladders["full"]
    gaps = {n: [rel(refined(table, S, full, l)[n][1]) for l in range(LEVELS)] for n in (11391, 11737)}
    for n, g in gaps.items():
        steps = [a - b for a, b in zip(g, g[1:])]
        assert all(s > 0 for s in steps) and all(0.45 < b / a < 0.55 for a, b in zip(steps, steps[1:]))
        floor = g[-1] - steps[-1]                       # geometric extrapolation (numerical)
        share = floor / g[0]
        print(n, "relative width by level", [f"{x:.5f}" for x in g], "extrapolated floor", f"{floor:.5f}",
              "share of level-0 width", f"{share:.3f}")
        lo_share, hi_share = {11391: (0.908, 0.910), 11737: (0.828, 0.829)}[n]
        assert lo_share < share < hi_share
    g = gaps[11391]
    assert 0.01017 < g[0] < 0.01019 and 0.00936 < g[-1] < 0.00937 and 0.920 < g[-1] / g[0] < 0.921
    g = gaps[11737]
    assert 0.00705 < g[0] < 0.00707 and 0.00599 < g[-1] < 0.00601 and 0.849 < g[-1] / g[0] < 0.850
    center = ladders["center"]
    cg = rel(refined(table, S, center, LEVELS - 1)[11391][1])
    assert 0.00597 < cg < 0.00598
    print("11391 width at level 3, eight other Yukawa inputs fixed at the center:", f"{cg:.5f}")

    # Width against the size of the eight-Yukawa box (one corner-anchored nested family).
    q = model.q
    def nested(frac):
        out = list(model.logs)
        for a in NONTOP:
            y = model.logs[a]
            lo = arb(y.lower())
            out[a] = lo.union(lo + (arb(y.upper()) - lo) * frac)
        return out
    scaling = {}
    for frac, (fa, fb) in ((q(1), ("2.20e-6", "2.22e-6")), (q("1/2"), ("1.68e-6", "1.70e-6")),
                           (q("1/4"), ("1.42e-6", "1.44e-6")), (q("1/16"), ("1.22e-6", "1.24e-6")),
                           (q(0), ("1.16e-6", "1.18e-6"))):
        p = vol["prepare"](model.gl, model.gh, nested(frac))
        face = p.faces[0][0].rad()
        lo_, hi_ = model.score(11391, p, 1024, 4)
        scaling[float(frac.mid())] = (hi_ - lo_) / hi_
        print("eight-Yukawa box fraction", frac.str(4), "first face radius", face.str(4),
              "11391 relative width at level 2", f"{scaling[float(frac.mid())]:.5f}")
        assert arb(fa) < face < arb(fb)
    w1, wh, w0 = scaling[1.0], scaling[0.5], scaling[0.0]
    assert 0.00947 < w1 < 0.00949 and 0.00778 < wh < 0.00780 and 0.00608 < w0 < 0.00610
    for lam, w in scaling.items():                  # the excess over the point box is linear in size
        assert abs(w - (w0 + lam * (w1 - w0))) < 1e-5, (lam, w)
    assert 0.00693 < scaling[0.25] < 0.00695 and 0.00629 < scaling[0.0625] < 0.00631

    # Corner probe (numerical): midpoints of V with the eight inputs fixed at four points.
    corner_mid = {}
    for tag in ("low", "high", "mixed", "center"):
        logs = list(model.logs)
        for j, a in enumerate(NONTOP):
            y = model.logs[a]
            pick = {"low": 0, "high": 1, "mixed": j % 2, "center": None}[tag]
            logs[a] = arb(y.mid()) if pick is None else arb(y.lower() if pick == 0 else y.upper())
        p = vol["prepare"](model.gl, model.gh, logs)
        corner_mid[tag] = {n: sum(model.score(n, p, 1024, 4)) / 2 for n in (11391, 11387, 11596, 11592, 11744, 11748)}
    for n in (11391, 11387, 11596, 11592, 11744, 11748):
        mids = [corner_mid[t][n] for t in corner_mid]
        spread = (max(mids) - min(mids)) / min(mids)
        assert spread < 5e-6, (n, spread)
    for (n, m), e in (((11391, 11387), -7.059e-3), ((11596, 11592), -4.658e-3), ((11744, 11748), 1.924e-3)):
        for t in corner_mid:
            d = (corner_mid[t][m] - corner_mid[t][n]) / corner_mid[t][n]
            assert abs(d - e) < 5e-7, (n, m, t, d)
    print("corner probe: midpoint spread of V below 5e-6 relative; pair midpoint differences",
          "-7.059e-3, -4.658e-3, 1.924e-3 at every corner")

    # Probes (first face only, full box): a narrower time band; one bisection of the gauge box.
    G = vol["prepare"].__globals__
    anchor, radius = G["anchor"], G["radius"]
    shift = q("-0.0183")
    try:
        G["radius"], G["anchor"] = q("1/128"), anchor + shift
        narrow = vol["prepare"](model.gl, model.gh, model.logs)
    finally:
        G["radius"], G["anchor"] = radius, anchor
    wide_face = model.preps["full"].faces[0][0].rad()
    narrow_face = narrow.faces[0][0].rad()
    print("face radius, band 1/8:", wide_face.str(4), "band 1/128:", narrow_face.str(4))
    assert arb("2.20e-6") < wide_face < arb("2.22e-6") and arb("2.12e-6") < narrow_face < arb("2.14e-6")
    splits = []
    for idx in itertools.product(range(2), repeat=3):
        a = [model.gl[i] + (model.gh[i] - model.gl[i]) * idx[i] / 2 for i in range(3)]
        b = [model.gl[i] + (model.gh[i] - model.gl[i]) * (idx[i] + 1) / 2 for i in range(3)]
        splits.append(vol["prepare"](a, b, model.center_logs))
    target, w = model.target(11391)
    lo_s = hi_s = arb(0)
    for p in splits:
        a, b = parts(target, p, arb(0), vol["radius"], 32)
        lo_s, hi_s = lo_s + a, hi_s + b
    a, b = parts(target, model.preps["center"], arb(0), vol["radius"], 32)
    one, eight = ((b - a) / b), ((hi_s - lo_s) / hi_s)
    print("11391 at the central top value: relative width with one gauge box", one.str(4),
          "with eight sub-boxes", eight.str(4))
    assert arb("0.0038") < one < arb("0.0039") and arb("0.0069") < eight < arb("0.0070")

    # 3. Classification at each level; previously decided labels never move.
    expect_pareto = {
        "full": [{}, {11722: "dominated"}, {11722: "dominated"}, {11722: "dominated"}],
        "center": [{10903: "dominated", 10972: "optimal", 11391: "optimal", 11722: "dominated"},
                   {10903: "dominated", 10972: "optimal", 11391: "optimal", 11722: "dominated",
                    11733: "dominated"},
                   {10903: "dominated", 10972: "optimal", 11391: "optimal", 11596: "optimal",
                    11722: "dominated", 11733: "dominated"},
                   {10903: "dominated", 10972: "optimal", 11391: "optimal", 11596: "optimal",
                    11722: "dominated", 11733: "dominated"}]}
    expect_lose = {"full": [[], [11273], [11273], [11273]],
                   "center": [[10918, 11015, 11172, 11273, 11322]] * 2
                   + [[10918, 11015, 11172, 11273, 11322, 11744]] * 2}
    for mode in ("full", "center"):
        for l in range(LEVELS):
            alt = refined(table, S, ladders[mode], l)
            lo, hi = R2["arrays"](alt, LABELS)
            st = dict(zip(LABELS, R2["pareto"](lo, hi, MAIN)))
            for n in LABELS:
                if st0[n] != "undecided":
                    assert st[n] == st0[n]
            decided = {n: st[n] for n in PARETO_UND if st[n] != "undecided"}
            assert decided == expect_pareto[mode][l], (mode, l, decided)
            exact = R2["fractions"](lo, hi)
            lin = {n: R2["selectability"](lo, hi, MAIN, row[n], exact)[0] for n in LINEAR_UND}
            assert [n for n in LINEAR_UND if lin[n] == "win"] == []
            assert [n for n in LINEAR_UND if lin[n] == "lose"] == expect_lose[mode][l], (mode, l)
            print(mode, "level", l, "Pareto decided", decided, "| linear: lose",
                  [n for n in LINEAR_UND if lin[n] == "lose"])

        # Full re-verification of the published counts at the finest level.
        alt = refined(table, S, ladders[mode], LEVELS - 1)
        lo, hi = R2["arrays"](alt, LABELS)
        st = dict(zip(LABELS, R2["pareto"](lo, hi, MAIN)))
        res = classify_linear(lo, hi, MAIN, LABELS, st)
        wins = [n for n, r in res.items() if r == "win"]
        assert wins == WINS_MAIN
        assert all(res.get(n, "lose") == "lose" for n in LABELS
                   if n not in LINEAR_UND and n not in WINS_MAIN)
        st5 = dict(zip(LABELS, R2["pareto"](lo, hi, FIVE)))
        wins5 = [n for n, r in classify_linear(lo, hi, FIVE, LABELS, st5).items() if r == "win"]
        assert len(wins5) == 34
        loa, hia = R2["arrays"](alt, OMEGA_A)
        sta = dict(zip(OMEGA_A, R2["pareto"](loa, hia, MAIN)))
        winsa = [n for n, r in classify_linear(loa, hia, MAIN, OMEGA_A, sta).items() if r == "win"]
        assert winsa == sorted(set(WINS_MAIN) - {10821, 10839, 10854} | {10875})
        assert minimax(lo, hi, LABELS, MAIN) == (11618, (175, 175))
        assert minimax(loa, hia, OMEGA_A, MAIN) == (11618, (163, 163))
        counts = {s: sum(v == s for v in st.values()) for s in ("dominated", "optimal", "undecided")}
        expect = {"full": {"dominated": 1169, "optimal": 102, "undecided": 7},
                  "center": {"dominated": 1171, "optimal": 105, "undecided": 2}}[mode]
        assert counts == expect
        print(mode, "finest level: Pareto", counts, "| linear wins 29, five-score wins 34, Omega_A wins 27,"
              " worst rank 175 / 163 unchanged")

    # 4. Overlaps that remain: per-level sequence and geometric limit (numerical).
    for mode in ("full", "center"):
        seqs = {}
        for n, m in PAIRS:
            seq = []
            for l in range(LEVELS):
                alt = refined(table, S, ladders[mode], l)
                (a, b), (c, d) = alt[n][1], alt[m][1]
                seq.append((min(b, d) - max(a, c)) / b)
            seqs[n, m] = seq
        for (n, m), seq in seqs.items():
            if seq[-1] > 0:
                steps = [a - b for a, b in zip(seq, seq[1:])]
                limit = seq[-1] - steps[-1] * steps[-1] / (steps[-2] - steps[-1])
                print(mode, n, m, "relative overlap", [f"{x:.2e}" for x in seq],
                      "geometric limit", f"{limit:.2e}")
                assert limit > 0.5 * seq[-1]
                assert all(0.49 < b / a < 0.51 for a, b in zip(steps, steps[1:]))
                e0, e3, el = EXPECT_OVERLAP[mode][n, m]
                for x, e in ((seq[0], e0), (seq[-1], e3), (limit, el)):
                    assert abs(x - e) < 0.005 * e, (mode, n, m, x, e)
        if mode == "full":
            assert [p for p, s in seqs.items() if s[-1] > 0] == [p for p in PAIRS if p != (11722, 11748)]
            s = seqs[11391, 11387]
            assert 0.00314 < s[0] < 0.00316 and 0.00233 < s[-1] < 0.00235
        else:
            assert [p for p, s in seqs.items() if s[-1] > 0] == [(11733, 11748), (11744, 11748)]
            s = seqs[11744, 11748]
            assert 0.00238 < s[0] < 0.00239 and 0.00133 < s[-1] < 0.00134

    # 5. Bisected eight-Yukawa box (256 sub-boxes): per-sub-box comparisons.
    assert sorted(subs) == list(range(256))
    for b in (0, 85, 170, 255):
        assert subbox_scores(model, b) == subs[b], b
    lvl3 = refined(table, S, ladders["full"], LEVELS - 1)

    def cut(b, n):
        (a, c), (x, y) = subs[b][n], lvl3[n][1]
        assert max(a, x) <= min(c, y)
        return max(a, x), min(c, y)

    m1 = min((cut(b, 10908)[0] - cut(b, 10903)[1]) / cut(b, 10903)[1] for b in subs)
    m2 = min((cut(b, 10972)[0] - cut(b, 10968)[1]) / cut(b, 10972)[1] for b in subs)
    print("sub-boxes: min relative margin V(10908)-V(10903)", f"{m1:.4e}", "V(10972)-V(10968)", f"{m2:.4e}")
    assert m1 > 0 and m2 > 0
    assert 0.01391 < m1 < 0.01393 and 0.000803 < m2 < 0.000804
    alt = refined(table, S, ladders["full"], LEVELS - 1)
    for n in SUB_LABELS:
        alt[n][1] = (min(cut(b, n)[0] for b in subs), max(cut(b, n)[1] for b in subs))
    lo, hi = R2["arrays"](alt, LABELS)
    st = dict(zip(LABELS, R2["pareto"](lo, hi, MAIN)))
    for n in LABELS:
        if st0[n] != "undecided":
            assert st[n] == st0[n]
    decided = {n: st[n] for n in PARETO_UND if st[n] != "undecided"}
    print("union over sub-boxes, Pareto decided:", decided)
    U = {n: alt[n][1] for n in SUB_LABELS}
    u1 = (U[10908][0] - U[10903][1]) / U[10903][1]
    u2 = (U[10972][0] - U[10968][1]) / U[10972][1]
    print("union margins: V(10908) over V(10903)", f"{u1:.7f}", "V(10972) over V(10968)", f"{u2:.7f}")
    assert 0.013855 < u1 < 0.013856 and 0.000783 < u2 < 0.000784
    assert 0 < u1 < m1 and 0 < u2 < m2
    for n in SUB_LABELS:                             # the union lies inside the level-3 interval
        assert lvl3[n][1][0] <= U[n][0] <= U[n][1] <= lvl3[n][1][1]
    # Test of the linear size law at k = 2 against the measured union overlaps (numerical).
    cen3 = refined(table, S, ladders["center"], LEVELS - 1)
    def rel_ov(t, n, m):
        (a, b), (c, d) = t[n][1], t[m][1]
        return (min(b, d) - max(a, c)) / b
    for (n, m), (pa, pb), (ma, mb) in (((10972, 10968), (-1.07e-3, -1.05e-3), (-7.85e-4, -7.83e-4)),
                                     ((10903, 10908), (-1.51e-2, -1.49e-2), (-1.386e-2, -1.385e-2))):
        pred = rel_ov(cen3, n, m) + (rel_ov(lvl3, n, m) - rel_ov(cen3, n, m)) / 2
        meas = rel_ov(alt, n, m)
        print("k = 2 test", n, m, "predicted", f"{pred:.4e}", "measured", f"{meas:.4e}",
              "optimism", f"{(meas - pred) / pred:.3f}")
        assert pa < pred < pb and ma < meas < mb and pred < meas < 0
    assert decided == {10903: "dominated", 10972: "optimal", 11722: "dominated"}
    counts = {k: sum(v == k for v in st.values()) for k in ("dominated", "optimal", "undecided")}
    assert counts == {"dominated": 1170, "optimal": 103, "undecided": 5}

    # Predicted bisection factor per axis for the rest, from the linear size law (numerical).
    cen = refined(table, S, ladders["center"], LEVELS - 1)
    need = {}
    for n, m in [(11391, 11387), (11596, 11592), (11733, 11726), (11744, 11722), (11744, 11748), (11733, 11748)]:
        ov = []
        for t in (lvl3, cen):
            (a, b), (c, d) = t[n][1], t[m][1]
            ov.append((min(b, d) - max(a, c)) / b)
        full_ov, point_ov = ov
        need[n, m] = math.floor((full_ov - point_ov) / -point_ov) + 1 if point_ov < 0 else None
        print("pair", n, m, "level-3 overlap", f"{full_ov:.3e}", "point limit", f"{point_ov:.3e}",
              "predicted factor per axis", need[n, m])
        expect_ov = {(11391, 11387): (2.340e-3, -1.066e-3), (11596, 11592): (2.303e-3, -1.590e-4),
                     (11733, 11726): (1.222e-3, -8.740e-4), (11744, 11722): (1.356e-3, -3.765e-4),
                     (11744, 11748): (2.723e-3, 1.335e-3), (11733, 11748): (3.752e-3, 2.026e-3)}[n, m]
        assert all(abs(x - e) < 0.005 * abs(e) for x, e in zip((full_ov, point_ov), expect_ov)), (n, m)
        if (n, m) == (11733, 11726):
            pred3 = point_ov + (full_ov - point_ov) / 3
            assert -1.76e-4 < pred3 < -1.74e-4
    assert need == {(11391, 11387): 4, (11596, 11592): 16, (11733, 11726): 3, (11744, 11722): 5,
                    (11744, 11748): None, (11733, 11748): None}

    # 10871: its upper bound stays positive while its sufficient window stays empty.
    for mode, cap in (("full", "3.43e-13"), ("center", "2.01e-13")):
        a, b = refined(table, S, ladders[mode], LEVELS - 1)[10871][1]
        assert a == 0.0 and 0 < b < float(cap)
        assert refined(table, S, ladders[mode], LEVELS - 1)[10854][1] == (0.0, 0.0)
    # Longest certified sufficient window of 10871 over 2048 top bins: still negative.
    target, _ = model.target(10871)
    pieces = 2048
    edges = [model.y0 + (model.y1 - model.y0) * j / pieces for j in range(pieces + 1)]
    for mode, (a, b) in (("full", ("-1.27e-6", "-1.26e-6")), ("center", ("-7.29e-7", "-7.28e-7"))):
        best = None
        for e0, e1 in zip(edges, edges[1:]):
            m = sufficient_margin(target, model.preps[mode], e0.union(e1) - model.logs[2].mid(), vol["radius"])
            best = m if best is None else best.max(m)
        print("10871", mode, "longest sufficient window", best.str(5))
        assert arb(a) < best < arb(b)
    print("total", round(time.time() - start), "s")


if __name__ == "__main__":
    args = sys.argv[1:]
    path = Path(args[args.index("--output") + 1]) if "--output" in args else DATA
    if "--regenerate" in args:
        regenerate(path)
        if "subboxes" not in json.loads(path.read_text()):
            print("no sub-box block in", path, "- run --subboxes 0 256 --output PATH, then verify")
            sys.exit(0)
    if "--subboxes" in args:
        i = args.index("--subboxes")
        generate_subboxes(path, int(args[i + 1]), int(args[i + 2]))
    else:
        verify(path)
