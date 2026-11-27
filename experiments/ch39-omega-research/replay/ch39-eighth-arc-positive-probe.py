"""第八壳最短弧四分之三处两种读出都能达到规范三式和固定能标正色。"""

from inspect import getclosurevars
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, ctx

SEEDS = (
    (
        "55.162111121792115",
        "56.077434945039712",
        "31.65946261659801",
        "69.169812458641516",
        "39.190363745416178",
        "77.773556357306376",
        "2.7781902884852103",
        "29.503657417330349",
        "28.999249962537689",
        "0.82719992028207867",
        "70.094482542775253",
        "77.79514028716018",
        "24.431945043575993",
        "67.565961934501431",
        "0.26728840283780625",
        "0.38068171855349892",
        "51.288527513920464",
        "14.315771947208523",
        "4.9672100498477079",
    ),
    (
        "54.939800217911305",
        "76.677456433191196",
        "67.652932098013167",
        "76.060363113952121",
        "0.27785755655752448",
        "75.99887719305643",
        "59.310934079710243",
        "0.93196284752615233",
        "6.175207118818955",
        "0.32920339015296124",
        "71.458951295571453",
        "1.2536561679824454",
        "8.9129793850988328",
        "64.255802283027549",
        "76.709512265045461",
        "0.42861530629175987",
        "51.756978284837878",
        "51.791048128818971",
        "26.368176144538641",
    ),
)
ACTIVE = ((15, 16, 11, 0, 18), (15, 16, 14, 0, 8))


def verify(root):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    arc = lib("eighth-shortest-arc-pivot-support")
    rows, single = arc["polynomials"](root)
    spectra = (single[1], arc["angle_weights"](rows[1], arb(3) / 4))
    num = arc["NUMBER"]
    _, source = lib("adjacent-ordered-charge-six-cutoffs")["verify"](
        root, refine_node=2_306_888_100_086
    )[2]
    source = [v.sqrt() for v in source]
    gauge = lib("all-field-signed-crossing")
    data = gauge["joint_system"](root, "1.176397e14", number=num, order=4, last=8)
    signs, masses = lib("rotation-portal-top-energy")["rotation_signs"](root, num)
    moments, heat, (field,) = lib("three-index-free-atoms")["broken_fields"](
        root, [num]
    )
    ctx.prec = 192
    neutral = lib("full-neutral-free-cutoff")
    species, _ = neutral["response_rows"](
        data,
        masses,
        signs,
        (arb("246.22") / arb("1.22089e19")) ** 2 / 2,
        arb(".42995401872"),
    )
    p0, p1 = data["p0"], data["p1"]
    probe = (arb(560000) / arb("1.22089e19")) ** 2
    pi = arb.pi()
    coeff = [
        (
            kind,
            z,
            ((q - 8 * c / arb(3)) / (48 * pi), (m - c) / (12 * pi), arb(c) / (40 * pi)),
        )
        for kind, z, (q, m, c) in species
    ]
    a1, a2, a3 = data["initial"]
    start = [(5 * a1 / 3 + a2 - 8 * a3 / 3) / 12, (a2 - a3) / 3, 0, (a3 - 2) / 10]
    weak = lib("native-rh-ray-pvm")["weak_series"](56)
    poly = getclosurevars(weak).nonlocals["coefficients"]
    means = lib("fourth-third-pivot-cyclic-five")["angular_means"]
    rule = gauge["angular_rule"](17, 45)
    ref = arb("-.074256577194663365545")

    def values(t):
        x = p0 * t.exp()
        sums = [arb(0)] * 3
        for kind, z, b in coeff:
            kernel = gauge["angular_exposure"](
                kind, z, x, p0, p1, (p1 / p0).log(), rule
            )
            sums[0] -= b[0] * kernel
            sums[1] -= b[1] * kernel
            if not b[2].is_zero():
                sums[2] -= b[2] * gauge["angular_exposure"](
                    kind, z, x, p0, probe, (probe / p0).log(), rule
                )
        d = weak(x / p0) if x < p0 * 10**7 else neutral["weak_column"](x, p0, p1)
        return sums[:2] + [100 * d, sums[2]]

    def slopes(t):
        x = p0 * t.exp()
        zero = means(p0 / x)
        end = [u - v for u, v in zip(zero, means(p1 / x), strict=True)]
        middle = [u - v for u, v in zip(zero, means(probe / x), strict=True)]
        sums = [arb(0)] * 3
        for kind, z, b in coeff:
            ex = (-z / x).exp()
            for i, ang in enumerate((end, end, middle)):
                sums[i] -= b[i] * ex * ang[kind]
        wd = (
            -sum(n * a * (p0 / x) ** n for n, a in enumerate(poly, 1))
            if x < p0 * 10**7
            else arb(0)
        )
        return sums[:2] + [100 * wd, sums[2]]

    def solve(k, ws, seed, axes):
        assert sum(ws).contains(1)
        tt = list(map(arb, seed))
        vals = [values(t) for t in tt]

        def residual():
            v = [
                b + sum(w * row[j] for w, row in zip(ws, vals, strict=True))
                for j, b in enumerate(start)
            ]
            return arb_mat(
                [
                    [u]
                    for u in [
                        *v,
                        sum(p0.sqrt() * (t / 2).exp() for t in tt) - sum(source),
                    ]
                ]
            )

        def jac(tlist):
            return arb_mat(
                [
                    [
                        *(ws[i] * v for v in slopes(tlist[i])),
                        p0.sqrt() * (tlist[i] / 2).exp() / 2,
                    ]
                    for i in axes
                ]
            ).transpose()

        for _ in range(2):
            shift = jac(tt).mid().solve(residual().mid()).entries()
            for i, d in zip(axes, shift, strict=True):
                tt[i] = (tt[i] - d).mid()
                vals[i] = values(tt[i])
        radius = arb("1e-5")
        box = [t + arb(0, radius.upper()) if i in axes else t for i, t in enumerate(tt)]
        pre = jac(tt).mid().inv().mid()
        assert not pre.det().contains(0)
        f = residual()
        noise = arb(0)
        for i in axes:
            x = p0 * box[i].exp()
            err = (
                arb("1e-19") + 2 * x / p1
                if x < p0 * 10**7
                else abs(neutral["weak_column"](x, p0, p1) - ref)
            )
            noise += abs(ws[i]) * err
        f[2, 0] += arb(0, (200 * noise).upper())
        diff = (pre * f).entries()
        error = (
            arb_mat([[int(i == j) for i in range(5)] for j in range(5)])
            - pre * jac(box)
        ).tolist()
        ratio = max(
            (abs(d) / radius + sum(map(abs, row))).upper()
            for d, row in zip(diff, error, strict=True)
        )
        gaps = lib("schur-horn-gravity-holdout")["ky_fan_gaps"](
            source, [(p0 * t.exp()).sqrt() for t in box]
        )
        thermal = [moments(field, heat, p0 * t.exp()) for t in box]
        gravity = [
            sum(w * v[i] for w, v in zip(ws, thermal, strict=True)) for i in (0, 1)
        ]
        print(
            "读出、五维自映比、主序最小余量、K、W",
            k,
            ratio.str(13),
            min(gaps).str(13),
            *[v.str(13) for v in gravity],
            flush=True,
        )
        limit, gap, kl, kh, wl, wh = (
            (".001322", "9.70e-20", ".06417", ".06418", "-.09670", "-.09668"),
            (".000555", "1.42e-19", ".07559", ".07560", "-.03485", "-.03483"),
        )[k]
        assert ratio < arb(limit) and min(gaps) > arb(gap)
        assert arb(kl) < gravity[0] < arb(kh) and arb(wl) < gravity[1] < arb(wh)

    for k, (ws, seed, axes) in enumerate(zip(spectra, SEEDS, ACTIVE, strict=True)):
        solve(k, ws, seed, axes)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
