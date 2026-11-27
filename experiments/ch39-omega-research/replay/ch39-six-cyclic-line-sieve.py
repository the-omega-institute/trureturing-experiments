"""严格三循环等变输运若能连续六壳满足三式，指定的低线选法仅余一种。"""

from bisect import bisect_left, bisect_right
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, ctx


def adj(m):
    return m.conjugate().transpose()


def squared(m):
    return sum(abs(z) ** 2 for z in m.entries())


def row(m, i):
    return acb_mat([[m[i, k] for k in range(m.ncols())]])


class Defects:
    def __init__(self, field):
        self.x = [arb(0)] + [arb(2) ** (arb(k) / 4) for k in range(-232, 1, 2)]
        self.parts = []
        for a in self.x:
            t = a * a
            vals = [arb(0) for _ in range(4)]
            if t > 0:
                for g, alpha, z in field:
                    if z.is_zero():
                        continue
                    u = z / t
                    if u > 140:
                        e2 = e3 = arb(0, arb("1e-60"))
                    else:
                        e2 = u.lower().expint(2).union(u.upper().expint(2))
                        e3 = u.lower().expint(3).union(u.upper().expint(3))
                    # K=-ck*a² - Σ alpha*t*(1-E2)/(2π).
                    # W= cw*a⁴ + Σ g*t²*(1/2-E3)/(4π).
                    vals[int(alpha > 0)] += abs(alpha) * t * (1 - e2) / (2 * arb.pi())
                    vals[2 + int(g < 0)] += (
                        abs(g) * t * t * (arb(1) / 2 - e3) / (4 * arb.pi())
                    )
            self.parts.append(vals)

    def interval(self, lo, hi):
        # 两端各向外取一个格点。热核亏损在相邻格点间单调。
        i = max(0, bisect_left(self.x, lo) - 1)
        j = min(len(self.x) - 1, bisect_right(self.x, hi))
        assert i == 0 or self.x[i] <= lo
        assert self.x[j] >= hi
        a = self.parts[i]
        b = self.parts[j]
        return [
            (a[k] - b[k + 1]).lower().union((b[k] - a[k + 1]).upper()) for k in (0, 2)
        ]


def thermal_defect(field, amp):
    # 用 1-exp(-u) <= min(1,u)，在 [1,max(1,t/z)] 与余段分开积分。
    t = amp**2
    dk = arb(0)
    dw = arb(0)
    for g, a, z in field:
        if z.is_zero():
            continue
        assert z > 0
        dk += abs(a) * z * (1 + max(arb(0), (t / z).log().upper()))
        dw += abs(g) * t * z
    return dk / (2 * arb.pi()), dw / (4 * arb.pi())


def prepare(levels, blocks, shells):
    ctx.prec = 256
    traces = list(map(squared, blocks))
    main = max((1, 2), key=lambda k: traces[k].lower())
    tiny = 3 - main
    assert traces[tiny] < arb("1e-9")
    model = []
    for frames, w, field in shells:
        assignment = []
        for i in range(19):
            active = [
                k for k, f in enumerate(frames) if squared(row(f, i)) > arb("0.9")
            ]
            assert len(active) == 1
            assignment.append(active[0])
        triv = [j for j in range(19) if assignment[j] == 0]
        big = [j for j in range(19) if assignment[j] == main]
        assert len(triv) == 7 and len(big) == 6
        assert all(w[i] < 0 for i in range(19) if assignment[i] != 0)
        assert (sum(r[1] for r in field) + arb(103) / 6).contains(0)
        assert arb(sum(r[0] for r in field) + 310).contains(0)
        ck = -sum(r[1] for r in field) / (2 * arb.pi())
        cw = -sum(r[0] for r in field) / (8 * arb.pi())
        # 分别将各字符的谱底和有质量热核误差计入充分宽的区间。
        errk, errw = thermal_defect(field, levels[0] + traces[main])
        assert errk * sum(abs(w[i]) for i in big) < arb(".012")
        assert errw * sum(abs(w[i]) for i in big) < arb(".002")
        tiny1 = levels[0] + traces[tiny]
        smallw = sum(abs(w[i]) for i in range(19) if assignment[i] == tiny)
        assert smallw * tiny1**2 * sum(abs(r[1]) for r in field) < arb("1e-12")
        assert smallw * tiny1**4 * sum(abs(r[0]) for r in field) < arb("1e-12")
        # 因为谱截止上界小于 1，二、四次幂的谱底平移逐行至多为 4*floor。
        assert levels[-1] < 1
        assert max(4 * levels[0] * sum(abs(v) for v in w) * v for v in (ck, cw)) < arb(
            "1e-10"
        )
        model.append(
            {
                "triv": triv,
                "big": big,
                "w": w,
                "zero": frames[0],
                "main": frames[main],
                "ck": ck,
                "cw": cw,
                "defects": Defects(field),
            }
        )
    return traces[0], traces[main], model


def bounds(model, levels, tzero, T, ids):
    B = acb_mat([row(shell["zero"], i).tolist()[0] for shell, i in zip(model, ids)])
    I = adj(B) * (B * adj(B)).inv()
    null = acb_mat([[int(i == j) for j in range(7)] for i in range(7)]) - I * B
    k = max(range(7), key=lambda j: null[j, j].real.lower())
    n = acb_mat([[null[j, k]] for j in range(7)])
    n = n / squared(n).sqrt()
    floor = levels[0]
    eps = (arb(3000) / arb("1.22089e19") - floor).sqrt()
    cap = eps * sum(sum(abs(I[k, j]) ** 2 for k in range(7)).sqrt() for j in range(6))
    assert cap**2 < tzero
    out = []
    for s in model:
        w = s["w"]
        K = arb(0)
        W = arb(0)
        ck = s["ck"]
        cw = s["cw"]
        for i in s["triv"]:
            f = row(s["zero"], i)
            align = abs((f * n)[0, 0])
            err = eps * sum(abs(z) for z in (f * I).entries())
            low = (
                floor
                + max(arb(0), (align * (tzero - cap**2).sqrt() - err).lower()) ** 2
            )
            high = floor + (align * tzero.sqrt() + err) ** 2
            z = low.lower().union(high.upper())
            # 热核亏损的上界在各场质量为零时消失，且随能标单调。
            ek, ew = s["defects"].interval(low.lower(), high.upper())
            K += w[i] * (-ck * z**2 + ek)
            W += w[i] * (cw * z**4 + ew)
        sk = (1 - K + arb(0, arb(".012001"))) / (ck * T**2)
        sw = (W - arb("1e-122") + arb(0, arb(".002001"))) / (cw * T**4)
        out.append((sk, sw))
    return out


def basic_options(model, targets):
    intervals = []
    for s, (u, v) in zip(model, targets):
        w = [-s["w"][i] for i in s["big"]]
        if (
            u.upper() < 0
            or v.upper() < 0
            or max(arb(0), u.lower()) ** 2 > sum(w) * v.upper()
        ):
            return False
        k = max(range(6), key=lambda i: w[i].lower())
        other = [a for i, a in enumerate(w) if i != k]
        constraints = []
        for p, target in ((2, u), (4, v)):
            a = sum(o ** (-arb(1) / (p - 1)) for o in other) ** (1 - p)
            constraints.append((p, target, a, max(o.upper() for o in other)))
        possible = []
        nodes = [(arb(0), arb(1), 0)]
        while nodes:
            l, r, depth = nodes.pop()
            if any(
                w[k] * l**p + a * (1 - r) ** p > target.upper()
                or w[k] * r**p + b * (1 - l) ** p < target.lower()
                for p, target, a, b in constraints
            ):
                continue
            if depth == 10:
                possible.append(l.union(r))
            else:
                m = (l + r) / 2
                nodes.extend(((l, m, depth + 1), (m, r, depth + 1)))
        if not possible:
            return False
        interval = possible[0].union(possible[-1])
        spectral = row(s["main"], s["big"][k])
        for prev, tab in intervals:
            overlap = abs((spectral * adj(prev))[0, 0]).lower()
            gap = max(arb(0), (1 - overlap**2).upper()).sqrt()
            if (
                interval.lower() > tab.upper() + gap
                or tab.lower() > interval.upper() + gap
            ):
                return False
        intervals.append((spectral, interval))
    return True


def experiment(root, fifth, columns, levels):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    api = lib("seven-cyclic-weak-horizon")
    candidates = api["CANDIDATES"][:6]
    graph = lib("three-ray-chiral-boundary")
    native = lib("native-rh-ray-pvm")
    gram = native["gram_cholesky"](root)
    frames, _ = lib("cyclic-two-characters-nine")["characters"](graph["geo"]["sites"])
    positions = range(2_306_888_100_092, 2_306_888_100_098)
    fields = lib("three-index-free-atoms")["broken_fields"](root, positions)[2]
    floor, heavy, top = levels
    weighted = (
        fifth
        * columns
        * acb_mat([[(heavy - floor).sqrt(), 0], [0, (top - floor).sqrt()]])
    )
    blocks = [adj(f) * weighted for f in frames]

    squared_levels, _, _, _, cover = lib("ninth-weak-two-line-threshold")["verify"](
        root
    )
    ctx.prec = 256
    assert floor.overlaps(squared_levels[0].sqrt())
    mz = arb("91.1876")
    planck = arb("1.22089e19")
    far = arb(".000020")
    near = arb(".020483945")
    cover((arb(3000) / mz) ** 2, far)
    cover((floor * planck / mz) ** 2, near)
    shells = []
    for number, ids, field in zip(positions, candidates, fields, strict=True):
        _, basis = lib("adjacent-pvm-transport")["spectrum"](graph, 6, number)
        weights = [
            z.real
            for z in native["readout"](basis, range(19), [14, 8, 7, 6], gram[3], 3)
        ]
        assert all(not w.contains(0) for w in weights) and all(
            weights[i] > 0 for i in ids
        )
        upper = arb("-.074256577194663365545") + sum(
            w * (far if i in ids else near if w > 0 else arb("-1.4e-11"))
            for i, w in enumerate(weights)
        )
        assert upper < 0
        projections = [adj(basis) * f for f in frames]
        assert all(squared(row(projections[0], i)) > arb(".9") for i in ids)
        shells.append((projections, weights, field))

    tzero, total, model = prepare(levels, blocks, shells)
    cover_api = lib("six-cyclic-diagonal-cover")
    geometry = cover_api["prepare"](model)
    survivors = []
    count = 0
    for ids in product(*candidates):
        count += 1
        targets = bounds(model, levels, tzero, total, ids)
        if basic_options(model, targets) and cover_api["test"](geometry, targets)[0]:
            survivors.append(ids)
    print(
        "六壳选线总数、排除数、剩余选线（一基）",
        count,
        count - len(survivors),
        [tuple(i + 1 for i in item) for item in survivors],
        flush=True,
    )
    assert count == 324 and survivors == [(15, 15, 0, 0, 0, 0)]


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    source = run_path(str(root / "ch39-fifth-full-positive-joint.py"))
    experiment(root, *source["verify"](root, return_operator="spectral"))
