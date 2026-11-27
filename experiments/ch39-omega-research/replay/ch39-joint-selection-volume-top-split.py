"""Top-Yukawa split enclosure of the volume score V: all undecided Pareto labels decided.

The joint comparison in ch39-five-criterion-pareto-selection and its refinement ladder
ch39-joint-selection-volume-ladder leave 5 Pareto labels (10871, 11391, 11596, 11733,
11744) and 22 linear-weight labels undecided, all through the width of the V enclosure
of ch39-volume-contenders.  This script encloses the same V with the same certificate
(`prepare` + the sufficient/necessary window integral `parts`) applied to eight
top-Yukawa sub-boxes, with three changes of computation, none of which changes the
enclosed quantity:

  1. the variational equation is integrated with one Picard refinement of the in-step
     derivative (the `refinements=1` option of ch39-volume-contenders.variations);
  2. the reference trajectory of the box centre uses 1048576 instead of 65536 steps;
  3. the Jacobian determinant is enclosed both directly and as det(pre m)/det(pre)
     and the two enclosures are intersected.

The other eight Yukawa coordinates and the gauge box stay whole, so every interval is
uniform in them, exactly as in the published comparison.  For each label the new V
interval is intersected with the level-3 ladder interval and with the published union
of the 256 bisected sub-boxes.

Each sub-box is padded outward by 2**-100; V is the mean over 8 x 128 equal log-top bins
with the density bounded on 16 time slabs.

Default run: loads the data file, recomputes one of the eight sub-box preparations and
the V scores of the labels that decide the Pareto questions bit for bit, and redoes every classification and
assertion from the data.  Writing modes (they overwrite the data file unless
`--output PATH` is given):
  --preps K               preparation of top sub-box K (0..7);
  --scores FIRST STOP     V scores of labels LABELS[FIRST:STOP] from the stored preparations;
  --full-check            recompute all preparations and all scores and compare.
"""
import json
import math
import sys
import time
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, ctx

root = Path(__file__).parent
DATA = root / "ch39-joint-selection-volume-top-split-data.json"
LAD = run_path(str(root / "ch39-joint-selection-volume-ladder.py"))
R2, vol = LAD["R2"], LAD["vol"]
gauge, q, ratios, running = vol["gauge"], vol["q"], vol["ratios"], vol["running"]
LABELS, MAIN, FIVE, OMEGA_A = R2["LABELS"], R2["MAIN"], R2["FIVE"], R2["OMEGA_A"]
T, CENTER_STEPS, BINS, SLABS = 8, 1048576, 128, 16
EPS = q(1) / 2**100            # outward padding of each top sub-box
PARETO_UND = [10871, 11391, 11596, 11733, 11744]
RIVALS = {10871: [10854], 11391: [11387], 11596: [11592], 11733: [11726, 11748],
          11744: [11722, 11748]}
SLICE = sorted(set(PARETO_UND) | {m for v in RIVALS.values() for m in v})
LINEAR_UND = [n for n in LAD["LINEAR_UND"] if n != 11273]     # 22 labels
EXPECT_WIDTH = {11391: (0.001109, 0.001111), 11744: (0.000771, 0.000772)}
# relative margins: (n, m) with n optimal uses (V(n)_lo - V(m)_hi) / V(n)_lo, with n dominated
# (V(m)_lo - V(n)_hi) / V(n)_hi; (11744, 11722) shows that 11722 does not dominate 11744.
EXPECT_MARGIN = {(11391, 11387): (0.00595, 0.00596), (11596, 11592): (0.00375, 0.00376),
                 (11733, 11726): (0.00398, 0.00399), (11744, 11748): (0.00115, 0.00116),
                 (11744, 11722): (0.00305, 0.00306)}
EXPECT_LINEAR = {"win": [11290, 11301, 11387, 11672],
                 "lose": [10918, 10946, 11015, 11172, 11322, 11355, 11381, 11398, 11409, 11744],
                 "undecided": [10968, 11069, 11193, 11391, 11499, 11521, 11527, 11596]}
EXPECT_SUMMARY = {"four": (33, 8), "five": (38, 8), "Omega_A": (33, 6)}


class Model(LAD["Model"]):
    """The ladder's model without its two whole-box preparations."""

    def __init__(self):
        ctx.prec = 112
        self.q = q
        self.gl, self.gh, self.logs = vol["input_bounds"]()
        fixed = gauge["fixed"]
        self.common = ((fixed["geo"]["BOOK_PLANCK"] / 1000)**2 * 3 * arb(5).sqrt()
                       / (4 * arb.pi()) * (-2 * vol["anchor"]).exp())
        self.offset = q(fixed["OFFSET"])
        self.weight = vol["base"]["spectral"]["weight"]
        self.y0 = 2 * (gauge["yukawa_centers"][2] * q(".99")).log()
        self.y1 = 2 * (gauge["yukawa_centers"][2] * q("1.01")).log()


class Bounds:
    """The six fields of FlowBounds that `parts` reads."""
    KEYS = ("top_sensitivity", "target_ratios", "transformed_rates", "rates")

    def __init__(self, faces, top_sensitivity, target_ratios, transformed_rates, rates, determinant):
        self.faces, self.determinant = faces, determinant
        self.top_sensitivity, self.target_ratios = top_sensitivity, target_ratios
        self.transformed_rates, self.rates = transformed_rates, rates


def prepare(gauge_lo, gauge_hi, logs, refinements=1, center_steps=CENTER_STEPS):
    """ch39-volume-contenders.prepare with changes 1-3 of the module docstring."""
    anchor = vol["anchor"]
    low = [a.union(b) for a, b in zip(gauge_lo, gauge_hi)]
    x_ref = [(a + b) / 2 for a, b in zip(gauge_lo, gauge_hi)]
    y_ref = [y.mid() for y in logs]
    jac, rates = vol["variations"](low, logs, refinements=refinements)
    m = arb_mat([[jac[i, j] for j in range(3)] for i in range(3)])
    pre = arb_mat([[m[i, j].mid() for j in range(3)] for i in range(3)]).inv()
    direct, conditioned = m.det(), (pre * m).det() / pre.det()
    assert direct.overlaps(conditioned)
    det = direct.intersection(conditioned)
    assert det > 0
    end, _, _ = gauge["evolve"](x_ref, y_ref, anchor, center_steps, diagonal_quarks=True)
    column = lambda values: arb_mat([[v] for v in values])
    normalized, slope, center = pre * column(ratios), pre * column(rates), pre * column(end)
    remainder = pre * jac
    for i in range(3):
        remainder[i, i] -= 1
    inputs = [v - x for v, x in zip(low, x_ref)] + [y - ref for y, ref in zip(logs, y_ref)]
    top = 5
    residual = [sum(remainder[i, j] * inputs[j] for j in range(12) if j != top) for i in range(3)]
    faces = [(center[i, 0] + residual[i] + (gauge_hi[i] - x_ref[i]),
              center[i, 0] + residual[i] + (gauge_lo[i] - x_ref[i])) for i in range(3)]
    return Bounds(faces, [remainder[i, top] for i in range(3)], [normalized[i, 0] for i in range(3)],
                  [slope[i, 0] for i in range(3)], rates, det)


def sub_edges(model, k):
    return (model.y0 + (model.y1 - model.y0) * k / T, model.y0 + (model.y1 - model.y0) * (k + 1) / T)


def prepare_sub(model, k):
    e0, e1 = sub_edges(model, k)
    logs = list(model.logs)
    logs[2] = (e0 - EPS).union(e1 + EPS)
    return logs[2].mid(), prepare(model.gl, model.gh, logs)


def ser(x):
    m, e = x.mid().man_exp()
    rm, re = x.rad().man_exp()
    return [int(m), int(e), int(rm), int(re)]


def des(s):
    """Ball with the stored midpoint and a radius at least the stored one."""
    m, e, rm, re = s
    return arb(arb(m) * arb(2)**e, arb(rm) * arb(2)**re)


def ser_bounds(ref, P):
    out = {"ref": ser(ref), "faces": [[ser(a), ser(b)] for a, b in P.faces],
           "determinant": ser(P.determinant)}
    for key in Bounds.KEYS:
        out[key] = [ser(x) for x in getattr(P, key)]
    return out


def des_bounds(d):
    P = Bounds([(des(a), des(b)) for a, b in d["faces"]], *[[des(x) for x in d[k]] for k in Bounds.KEYS],
               des(d["determinant"]))
    return des(d["ref"]), P


def score(model, preps, n):
    """V interval of label n: T sub-boxes x BINS equal log-top bins, SLABS time slabs."""
    target, w = model.target(n)
    a, b = arb(0), arb(0)
    for k, (ref, P) in enumerate(preps):
        e0, e1 = sub_edges(model, k)
        edges = [e0 + (e1 - e0) * j / BINS for j in range(BINS + 1)]
        for x0, x1 in zip(edges, edges[1:]):
            lo, up = LAD["parts"](target, P, x0.union(x1) - ref, vol["radius"], SLABS)
            a, b = a + lo, b + up
    lo, up = (w * a / (T * BINS)).lower(), (w * b / (T * BINS)).upper()
    return (math.nextafter(float(lo), -math.inf) if lo > 0 else 0.0,
            math.nextafter(float(up), math.inf) if up > 0 else 0.0)


def fbits(x):
    return float.hex(x)


def write_preps(path, k):
    model = Model()
    raw = json.loads(path.read_text()) if path.exists() else {}
    block = raw.setdefault("preps", [None] * T)
    ref, P = prepare_sub(model, k)
    d = ser_bounds(ref, P)
    ref2, P2 = des_bounds(d)                  # decoding only widens a radius
    assert ref2.contains(ref) and P2.determinant.contains(P.determinant)
    assert all(b.contains(a) for x, y in zip(P.faces, P2.faces) for a, b in zip(x, y))
    assert all(b.contains(a) for key in Bounds.KEYS for a, b in zip(getattr(P, key), getattr(P2, key)))
    block[k] = d
    path.write_text(json.dumps(raw, sort_keys=True))
    print("sub-box", k, "prepared", flush=True)


def write_scores(path, first, stop):
    model = Model()
    raw = json.loads(path.read_text())
    preps = [des_bounds(d) for d in raw["preps"]]
    block = raw.setdefault("scores", {})
    for n in LABELS[first:stop]:
        block[str(n)] = [fbits(x) for x in score(model, preps, n)]
    path.write_text(json.dumps(raw, sort_keys=True))
    print("wrote scores", LABELS[first], "to", LABELS[stop - 1], flush=True)


def baseline():
    """Level-3 ladder intervals with the published sub-box unions (1170 / 103 / 5)."""
    table, S, ladders, subs = LAD["load"]()
    alt = LAD["refined"](table, S, ladders["full"], LAD["LEVELS"] - 1)
    for n in LAD["SUB_LABELS"]:
        lv = alt[n][1]
        cuts = [(max(subs[b][n][0], lv[0]), min(subs[b][n][1], lv[1])) for b in subs]
        alt[n][1] = (min(c[0] for c in cuts), max(c[1] for c in cuts))
    return table, alt


def verify(path=DATA, full=False):
    start = time.time()
    raw = json.loads(path.read_text())
    assert len(raw["preps"]) == T and sorted(map(int, raw["scores"])) == LABELS
    model = Model()
    preps = [des_bounds(d) for d in raw["preps"]]
    new = {int(n): tuple(float.fromhex(x) for x in v) for n, v in raw["scores"].items()}

    # 1. Bit-exact recomputation: one preparation (all eight with --full-check) and the
    #    scores of the labels that decide the Pareto questions (all with --full-check).
    for k in (range(T) if full else [1]):
        ref, P = prepare_sub(model, k)
        assert ser_bounds(ref, P) == raw["preps"][k], k
    print("preparations recomputed in", round(time.time() - start), "s", flush=True)
    check = LABELS if full else sorted(set(SLICE) | set(LINEAR_UND))
    for n in check:
        assert score(model, preps, n) == new[n], n
    print("scores of", len(check), "labels recomputed, total", round(time.time() - start), "s", flush=True)

    # 2. Where the width came from: on the whole box (no top split, 65536 centre steps,
    #    1024 bins, 4 slabs) one Picard refinement alone cuts the width of 11391 from the
    #    published level-2 value 0.00948 to 0.00547 and the first face radius by half.
    original = vol["prepare"](model.gl, model.gh, model.logs)
    refined1 = prepare(model.gl, model.gh, model.logs, refinements=1, center_steps=65536)
    w = {}
    for tag, P in (("original", original), ("refinement", refined1)):
        a, b = model.score(11391, P, 1024, 4)
        w[tag] = (b - a) / b
        print(tag, "first face radius", P.faces[0][0].rad().str(3), "11391 relative width", f"{w[tag]:.5f}")
    assert 0.00947 < w["original"] < 0.00949 and 0.00546 < w["refinement"] < 0.00548
    assert arb("2.20e-6") < original.faces[0][0].rad() < arb("2.22e-6")
    assert arb("1.07e-6") < refined1.faces[0][0].rad() < arb("1.09e-6")

    # 3. Width budget per sub-box: face radii, top sensitivity and determinant per sub-box.
    for k, (ref, P) in enumerate(preps):
        faces = max(x.rad() for pair in P.faces for x in pair)
        sens = max(x.rad() for x in P.top_sensitivity)
        if k in (0, T - 1):
            print("sub-box", k, "largest face radius", faces.str(3), "largest top-sensitivity radius",
                  sens.str(3), "determinant radius", P.determinant.rad().str(3))
        assert faces < arb("4.0e-7") and sens < arb("3.8e-5") and P.determinant.rad() < arb("1.4e-4")

    # 4. Intersect with the published intervals and reclassify (Pareto).
    table, alt = baseline()
    lo0, hi0 = R2["arrays"](alt, LABELS)
    st0 = dict(zip(LABELS, R2["pareto"](lo0, hi0, MAIN)))
    assert [n for n in LABELS if st0[n] == "undecided"] == PARETO_UND
    assert sum(v == "dominated" for v in st0.values()) == 1170
    row = {n: k for k, n in enumerate(LABELS)}
    for n, rv in RIVALS.items():
        k = row[n]
        poss = [LABELS[m] for m in range(len(LABELS)) if m != k and all(hi0[m, :4] >= lo0[k, :4])]
        assert poss == rv
        # every possible dominator is already certainly better in W, -eta and -E
        assert all(all(lo0[row[m], [0, 2, 3]] > hi0[k, [0, 2, 3]]) for m in rv)
    for n, (a, b) in new.items():
        c, d = alt[n][1]
        assert max(a, c) <= min(b, d), n
        alt[n][1] = (max(a, c), min(b, d))
    assert all(alt[n][1] == new[n] for n in SLICE)          # the new interval is the tighter one
    rel = lambda n: (alt[n][1][1] - alt[n][1][0]) / alt[n][1][1]
    for n, (lo_, hi_) in EXPECT_WIDTH.items():
        assert lo_ < rel(n) < hi_, (n, rel(n))
    print("relative width of V:", {n: f"{rel(n):.6f}" for n in EXPECT_WIDTH})

    lo, hi = R2["arrays"](alt, LABELS)
    st = dict(zip(LABELS, R2["pareto"](lo, hi, MAIN)))
    for n in LABELS:
        if st0[n] != "undecided":
            assert st[n] == st0[n]
    decided = {n: st[n] for n in PARETO_UND}
    print("Pareto status of the five:", decided)
    assert decided == {10871: "optimal", 11391: "optimal", 11596: "optimal",
                       11733: "dominated", 11744: "dominated"}
    counts = {s: sum(v == s for v in st.values()) for s in ("dominated", "optimal", "undecided")}
    assert counts == {"dominated": 1172, "optimal": 106, "undecided": 0}
    V = lambda n: alt[n][1]
    margins = {}
    for n, m in ((11391, 11387), (11596, 11592)):        # optimal: the rival certainly has smaller V
        margins[n, m] = (V(n)[0] - V(m)[1]) / V(n)[0]
    for n, m in ((11733, 11726), (11744, 11748)):        # the rival certainly has larger V
        margins[n, m] = (V(m)[0] - V(n)[1]) / V(n)[1]
    margins[11744, 11722] = (V(11744)[0] - V(11722)[1]) / V(11744)[0]
    for p, (a, b) in EXPECT_MARGIN.items():
        print("relative margin", p, f"{margins[p]:.6f}")
        assert a < margins[p] < b, (p, margins[p])
    assert V(11748)[1] < V(11733)[0]                     # 11748 certainly has smaller V: no domination
    assert V(10854) == (0.0, 0.0)
    assert 1.350e-14 < V(10871)[0] < 1.351e-14 and 4.42e-14 < V(10871)[1] < 4.43e-14
    print("10871: V in", V(10871), "and V(10854) = 0")

    # 5. 10871: top bins whose sufficient window is certified non-empty.
    target, _ = model.target(10871)
    positive, best = [], None
    for k, (ref, P) in enumerate(preps):
        e0, e1 = sub_edges(model, k)
        edges = [e0 + (e1 - e0) * j / BINS for j in range(BINS + 1)]
        for j, (x0, x1) in enumerate(zip(edges, edges[1:])):
            m = LAD["sufficient_margin"](target, P, x0.union(x1) - ref, vol["radius"])
            if m > 0:
                positive.append(k * BINS + j)
            best = m if best is None else best.max(m)
    print("10871: top bins with a certified non-empty sufficient window", positive[0], "..", positive[-1],
          "longest window", best.str(5))
    assert positive == list(range(104)) and arb("3.705e-7") < best < arb("3.706e-7")

    # 6. Linear-weight selectability of the 22 undecided labels; nothing decided moves.
    exact = R2["fractions"](lo, hi)
    lin = {}
    for n in LABELS:
        if st[n] == "dominated":
            continue
        lin[n] = R2["selectability"](lo, hi, MAIN, row[n], exact)[0]
    wins = [n for n, r in lin.items() if r == "win"]
    assert set(LAD["WINS_MAIN"]) <= set(wins)
    for n in LABELS:                                    # previously certified non-selectable stay so
        if n not in LAD["WINS_MAIN"] and n not in LAD["LINEAR_UND"]:
            assert lin.get(n, "lose") == "lose", n
    status = {n: ("dominated" if st[n] == "dominated" else lin[n]) for n in LINEAR_UND}
    print("linear status of the 22:", status)
    assert {n for n, r in status.items() if r == "win"} == set(EXPECT_LINEAR["win"])
    assert {n for n, r in status.items() if r in ("lose", "dominated")} == set(EXPECT_LINEAR["lose"])
    assert [n for n, r in status.items() if r == "undecided"] == EXPECT_LINEAR["undecided"]
    assert [n for n in LABELS if lin.get(n) == "undecided"] == EXPECT_LINEAR["undecided"]
    assert sorted(wins) == sorted(set(LAD["WINS_MAIN"]) | set(EXPECT_LINEAR["win"]))
    st5 = dict(zip(LABELS, R2["pareto"](lo, hi, FIVE)))
    lin5 = {n: R2["selectability"](lo, hi, FIVE, row[n], exact)[0] for n in LABELS if st5[n] != "dominated"}
    loa, hia = R2["arrays"](alt, OMEGA_A)
    sta = dict(zip(OMEGA_A, R2["pareto"](loa, hia, MAIN)))
    exa = R2["fractions"](loa, hia)
    rowa = {n: k for k, n in enumerate(OMEGA_A)}
    lina = {n: R2["selectability"](loa, hia, MAIN, rowa[n], exa)[0] for n in OMEGA_A if sta[n] != "dominated"}
    summary = {"four": (len(wins), sum(r == "undecided" for r in lin.values())),
               "five": (sum(r == "win" for r in lin5.values()), sum(r == "undecided" for r in lin5.values())),
               "Omega_A": (sum(r == "win" for r in lina.values()), sum(r == "undecided" for r in lina.values()))}
    print("linear winners / undecided:", summary)
    assert summary == EXPECT_SUMMARY
    old_a = set(LAD["WINS_MAIN"]) - {10821, 10839, 10854} | {10875}
    assert old_a <= {n for n, r in lina.items() if r == "win"}
    # Five scores on Omega and four scores on Omega_A: recompute the classification on the
    # published intervals and check that every win stays a win and every lose stays a lose.
    for cols, labs, (lb, hb), new_lin, new_st in ((FIVE, LABELS, (lo0, hi0), lin5, st5),
                                                  (MAIN, OMEGA_A, R2["arrays"](baseline()[1], OMEGA_A), lina, sta)):
        rr = {n: k for k, n in enumerate(labs)}
        stb = dict(zip(labs, R2["pareto"](lb, hb, cols)))
        exb = R2["fractions"](lb, hb)
        old = {n: ("lose" if stb[n] == "dominated" else R2["selectability"](lb, hb, cols, rr[n], exb)[0])
               for n in labs}
        now = {n: ("lose" if new_st[n] == "dominated" else new_lin[n]) for n in labs}
        assert all(now[n] == old[n] for n in labs if old[n] != "undecided")
        print("old wins", sum(r == "win" for r in old.values()), "kept; old non-selectable kept",
              sum(r == "lose" for r in old.values()))
        if cols == FIVE:
            assert sum(r == "win" for r in old.values()) == 34
        else:
            assert {n for n, r in old.items() if r == "win"} == old_a
    assert LAD["minimax"](lo, hi, LABELS, MAIN) == (11618, (175, 175))
    assert LAD["minimax"](loa, hia, OMEGA_A, MAIN) == (11618, (163, 163))
    # 7. Numerical probe only: shrink every V interval about its midpoint by a factor and
    #    rerun the linear test on the 8 undecided labels.  Not a certificate.
    probe = {}
    for factor in (0.0, 0.5, 0.75):
        l2, h2 = lo.copy(), hi.copy()
        mid, half = (lo[:, 1] + hi[:, 1]) / 2, (hi[:, 1] - lo[:, 1]) / 2 * factor
        l2[:, 1], h2[:, 1] = mid - half, mid + half
        ex2 = R2["fractions"](l2, h2)
        probe[factor] = {n: R2["selectability"](l2, h2, MAIN, row[n], ex2)[0] for n in EXPECT_LINEAR["undecided"]}
        print("probe, V half-widths x", factor, probe[factor])
    assert probe[0.0] == {10968: "win", 11069: "win", 11193: "win", 11391: "win", 11499: "win",
                          11521: "win", 11527: "win", 11596: "lose"}
    und = {n: "undecided" for n in EXPECT_LINEAR["undecided"]}
    assert probe[0.5] == {**und, 11193: "win", 11391: "win", 11499: "win", 11596: "lose"}
    assert probe[0.75] == {**und, 11391: "win", 11596: "lose"}
    print("worst rank: 11618 unique, 175 on Omega and 163 on Omega_A; total", round(time.time() - start), "s")


if __name__ == "__main__":
    args = sys.argv[1:]
    path = Path(args[args.index("--output") + 1]) if "--output" in args else DATA
    if "--preps" in args:
        write_preps(path, int(args[args.index("--preps") + 1]))
    elif "--scores" in args:
        i = args.index("--scores")
        write_scores(path, int(args[i + 1]), int(args[i + 2]))
    else:
        verify(path, full="--full-check" in args)
