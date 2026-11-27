"""重放预登记 9b57c857 的第六壳味留出，生成无选相的有理扰动证书输入。"""
from argparse import ArgumentParser
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, ctx


NUMBER = 2_306_888_100_091
SCALE = 10**14
TABLE = "ch39-sixth-quark-distance-cells.tsv"


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    w, dw, sectors, _, _ = lib("fifth-third-pivot-joint")["verify"](root, return_transport=True)
    assert dw < arb("2.81e-9")
    ctx.prec = 384
    shell = lib("three-ray-chiral-boundary")
    spec = lib("adjacent-pvm-transport")["spectrum"]
    adj = lib("rank-two-generator-up-quark")["adj"]
    windows = lib("three-shell-outer-holdout")["WINDOWS"]
    sites = (shell["soft"]["pole"], *shell["boundary_sites"]())
    assert sites == (14, 8, 7, 6)
    rows = [("scale", SCALE), ("shell", NUMBER), ("sites", *sites), ("dw_upper", 281000, SCALE)]

    def floor(v):
        return int((v*SCALE).lower().floor().unique_fmpz())

    def ceil(v):
        return int((v*SCALE).upper().ceil().unique_fmpz())

    def grid(v):
        k = int((v.mid()*SCALE+arb(1)/2).floor().unique_fmpz())
        assert abs(v*SCALE-k) < 1
        return k

    for iq, (charge, power, op) in enumerate(((-4, 2, adj(w**4)), (2, 1, w**2))):
        src, _, train, _, _ = sectors[iq]
        dst = spec(shell, charge, NUMBER)[1]
        overlap = adj(dst)*op*src
        probs = [[abs(overlap[i, j])**2 for j in range(19)] for i in range(19)]
        assert all(sum(p).contains(1) for p in probs)
        weights = [sum(a*b for a, b in zip(p, train, strict=True)) for p in probs]
        four = acb_mat([[dst[i, j] for j in range(19)] for i in sites])
        raw = lib("common-unitary-quark-right-trace")["weighted"](four, weights)
        # 对两次 W 操作使用同一算子误差；上游 Z 的整个盒已由 src 包围。
        error = 2*power*dw*(max(v.upper() for v in train)-min(v.lower() for v in train))
        levels = sorted((z.real for z in raw.eig(algorithm="rump")), key=lambda z: z.mid())
        levels = [v+arb(0, error.upper()) for v in levels]
        assert all(a < b for a, b in zip(levels, levels[1:]))
        status = []
        for j, (low, high) in enumerate(windows[iq]):
            a, b = arb(low), arb(high)
            if levels[j] > b or levels[j+1] < a:
                status.append("排除")
            elif levels[j] < a and levels[j+1] > b:
                status.append("覆盖全窗")
            else:
                status.append("未决")
        bands = (
            ((".00007757", ".00007763"), (".00375553", ".00375567"),
             (".51618076", ".51618084"), (".92965613", ".92965627")),
            ((".000054527", ".000054533"), (".00057938", ".00057940"),
             (".00225657", ".00225659"), (".09478857", ".09478859")),
        )
        assert all(arb(a) < v < arb(b) for v, (a, b) in zip(levels, bands[iq], strict=True))
        assert status == ["排除", "排除", "覆盖全窗"]
        print("荷、四站本征值、单窗判定", charge, *[v.str(13) for v in levels], status, flush=True)
        rows.append(("train", charge, *map(floor, train)))
        for i, prob in enumerate(probs):
            rows.append(("prob", charge, i, *map(ceil, prob)))
        for j in range(19):
            cells = []
            for i in range(4):
                for k in range(i, 4):
                    v = four[i, j]*four[k, j].conjugate()
                    cells.extend((grid(v.real), grid(v.imag)))
            rows.append(("projector", charge, j, *cells))
    print("训练盒中 W 至中心的算子误差上界", dw.str(12), flush=True)
    return "".join("\t".join(map(str, row))+"\n" for row in rows)


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="重建有限有理输入表")
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    data = verify(root)
    if args.write:
        (root/TABLE).write_text(data)
    else:
        assert (root/TABLE).read_text() == data
