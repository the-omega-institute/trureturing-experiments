"""固定第五壳源谱，跨字符的一次酉输运在第七至十二壳共同实现十八式。"""

from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat, ctx

LINES = (15, 15, 0, 0, 0, 0)
ACTIVE = (2, 5, 9, 12, 17, 24, 28, 33, 39, 41, 47, 48)
ENERGIES = (
    "954.6169703989897",
    "895.5633453452054",
    "544.705683198702",
    "544.755583187606",
    "533.1938263366604",
    "516.6903911117098",
)
PARAMETERS = (
    "1.625",
    "-0.125",
    "2.8936073969540357",
    "-3.875",
    "-1.625",
    "1.1406038926069186",
    "2.125",
    "-1.125",
    "1.125",
    "-2.0763812281846183",
    "2.5",
    "-0.25",
    "0.13750738286813155",
    "0.75",
    "-0.5",
    "0.5",
    "-2.0",
    "-0.9846155425712885",
    "0.375",
    "-0.5",
    "-1.0",
    "-0.125",
    "1.125",
    "0.625",
    "2.478428562359693",
    "1.75",
    "2.625",
    "-2.0",
    "-0.2277570371104276",
    "2.25",
    "1.125",
    "-0.5",
    "0.875",
    "-0.6789753340020642",
    "-0.25",
    "-0.625",
    "0.875",
    "2.5",
    "1.5",
    "0.5519513135317411",
    "-1.0",
    "-1.4847576987464315",
    "3.25",
    "0.5",
    "2.625",
    "0.125",
    "0.75",
    "-2.442684815545423",
    "-2.400132242260197",
    "0.75",
)


def col(m, i):
    return acb_mat([[m[k, i]] for k in range(m.nrows())])


def join(a, b):
    return acb_mat([i + j for i, j in zip(a.tolist(), b.tolist(), strict=True)])


def certify(root, levels, shells, source_trace):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    common = lib("common-three-shell-nine-endpoints")
    adj, inner, normalize = common["adj"], common["inner"], common["normalize"]
    sq = lib("six-cyclic-line-sieve")["squared"]
    gravity = lib("shortest-common-complex-gate")["gravity_system"]
    heat = lib("chiral-minimal-neutral-gravity")["expint"]
    ctx.prec = 256
    floor, first, second = levels
    gain = second - floor
    g2 = first - floor
    lines = LINES
    active = ACTIVE
    B = acb_mat([[f[0][i, k] for k in range(7)] for (f, w, z), i in zip(shells, lines)])
    J = adj(B) * (B * adj(B)).inv()
    F = acb_mat([[int(i == j) for j in range(7)] for i in range(7)]) - J * B
    n = col(F, 4)
    n = n / sq(n).sqrt()
    mp = arb("1.22089e19")
    mz = arb("91.1876")
    base = list(map(arb, PARAMETERS))
    C = [base[i] for i in active]
    Q = list(map(arb, ENERGIES))
    qwidth = arb("2e-8")
    qbox = [q + arb(0, qwidth.upper()) for q in Q]
    frames = [(f[0], join(f[1], f[2])) for f, _, _ in shells]

    def frame(x, qs):
        v = list(base)
        for i, t in zip(active, x):
            v[i] = t
        raw = acb_mat(
            [
                [acb(v[2 * i + j], v[24 + 2 * i + j]) for j in range(2)]
                for i in range(12)
            ]
        )
        a, b = col(raw, 0), col(raw, 1)
        ds = []
        for k in active:
            d = acb_mat(12, 2)
            if k < 48:
                im, t = divmod(k, 24)
                r, j = divmod(t, 2)
                d[r, j] = acb(1 - im, im)
            ds.append(d)
        u, du = normalize(a, [col(d, 0) for d in ds])
        s = inner(u, b)
        v2, dv = normalize(
            b - u * s,
            [
                col(d, 1) - z * s - u * (inner(z, b) + inner(u, col(d, 1)))
                for d, z in zip(ds, du)
            ],
        )
        y = acb_mat([[((q / mp - floor).sqrt())] for q in qs])
        lift = J * y
        m = sq(lift)
        g = (gain - m).sqrt()
        theta = v[48]
        phi = v[49]
        z = acb(0, phi).exp()
        a0 = g * theta.sin() * z
        a1 = g * theta.cos()
        z2 = g2.sqrt()
        c = join(lift + n * a0, acb_mat(7, 1))
        d = join(u * a1, v2 * z2)
        d0 = [
            join(
                n * (a1 * z if k == 48 else a0 * acb(0, 1) if k == 49 else 0),
                acb_mat(7, 1),
            )
            for k in active
        ]
        d1 = [
            join(
                uu * a1 - (u * (g * theta.sin()) if k == 48 else acb_mat(12, 1)),
                vv * z2,
            )
            for k, uu, vv in zip(active, du, dv)
        ]
        return c, d, d0, d1

    def systems(x, qs, values, read=False):
        c, d, dc, dd = frame(x, qs)
        R = []
        D = []
        pr = []
        for (f, w, field), (zero, merged) in zip(shells, frames):
            U = zero * c + merged * d
            if read:
                pr.append(U)
                continue
            derivatives = [zero * a + merged * b for a, b in zip(dc, dd)]
            a, j = gravity(
                U, derivatives, (floor, floor + 1, floor + 1), w, heat, field, values
            )
            R.extend(a.tolist())
            D.extend(j.tolist())
        return pr if read else (arb_mat(R), arb_mat(D))

    for _ in range(3):
        rr, j = systems(C, Q, True)
        step = j.mid().solve(rr.mid()).entries()
        C = [(a - b).mid() for a, b in zip(C, step)]
    f, j = systems(C, qbox, True)
    pre = j.mid().inv().mid()
    off = [v.abs_upper() for v in (pre * f).entries()]

    widths = [2 * e + arb("3e-9") for e in off]
    box = [a + arb(0, w.upper()) for a, w in zip(C, widths)]
    _, wide = systems(box, qbox, False)
    defect = arb_mat([[int(i == j) for j in range(12)] for i in range(12)]) - pre * wide
    ratio = max(
        (o + sum(abs(c).upper() * w for c, w in zip(row, widths))) / width
        for o, row, width in zip(off, defect.tolist(), widths)
    )
    print("十二维自映比与最大盒半宽", ratio.str(10), max(widths).str(12), flush=True)
    assert ratio < arb(".051") and not pre.det().contains(0)
    ref = arb("-.074256577194663365546").union(arb("-.074256577194663365544"))
    weak = lib("native-rh-ray-pvm")["weak_series"](56)
    large = lib("full-neutral-free-cutoff")["weak_column"]
    for (_, w, _), U, i, q in zip(
        shells, systems(box, qbox, False, True), lines, Q, strict=True
    ):
        amin = min(
            (floor + sum(abs(z) ** 2 for z in row)).lower()
            for j, row in enumerate(U.tolist())
            if j != i
        )
        p0 = (mz / mp) ** 2
        p1 = (arb("1.176397e14") / mp) ** 2
        assert amin**2 > p0
        tail = (large(amin**2, p0, p1) - ref).abs_upper() * sum(
            abs(v).upper() for j, v in enumerate(w) if j != i
        )
        F = [
            ref + w[i] * (weak(((q + sign * qwidth) / mz) ** 2) - ref) + arb(0, tail)
            for sign in (-1, 1)
        ]
        print(
            "其他线的截止下界与两个弱差面",
            amin.str(9),
            *[a.str(9) for a in F],
            flush=True,
        )
        assert F[0] > 0 > F[1]

    blocks = frame(box, qbox)
    image_trace = sq(blocks[0])
    assert source_trace > arb(".492") and image_trace < arb(".345")
    print(
        "平凡字符的源、像加权迹", source_trace.str(12), image_trace.str(12), flush=True
    )
    return blocks[:2]


def experiment(root, fifth, columns, levels):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    graph = lib("three-ray-chiral-boundary")
    native = lib("native-rh-ray-pvm")
    adj = lib("common-three-shell-nine-endpoints")["adj"]
    frames, _ = lib("cyclic-two-characters-nine")["characters"](graph["geo"]["sites"])
    floor, first, second = levels
    source = (
        fifth
        * columns
        * acb_mat([[(first - floor).sqrt(), 0], [0, (second - floor).sqrt()]])
    )
    source_trace = lib("six-cyclic-line-sieve")["squared"](adj(frames[0]) * source)
    gram = native["gram_cholesky"](root)
    positions = range(2_306_888_100_092, 2_306_888_100_098)
    fields = lib("three-index-free-atoms")["broken_fields"](root, positions)[2]
    shells = []
    for number, field in zip(positions, fields, strict=True):
        _, basis = lib("adjacent-pvm-transport")["spectrum"](graph, 6, number)
        trivial_line = (adj(basis) * frames[0]).tolist()[LINES[len(shells)]]
        assert sum(abs(z) ** 2 for z in trivial_line) > arb(".9")
        w = [
            v.real
            for v in native["readout"](basis, range(19), [14, 8, 7, 6], gram[3], 3)
        ]
        shells.append(([adj(basis) * f for f in frames], w, field))
    c, d = certify(root, levels, shells, source_trace)
    return frames[0] * c + join(frames[1], frames[2]) * d


if __name__ == "__main__":
    root = Path(__file__).resolve().parent
    src = run_path(str(root / "ch39-fifth-full-positive-joint.py"))
    experiment(root, *src["verify"](root, return_operator="spectral"))
