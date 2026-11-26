"""由三壳训练谱振幅按相邻完整投影递推，在两侧各检验两枚实际外壳。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, ctx


CENTER = 2_306_888_100_085
CHARGES = (-4, 2)


def training(root, sign, source):
    if sign < 0:
        geometry = (".0182", ".113", ("-.21668", "-.53516", ".71837", "-.38807"))
        return source(root, ((CENTER-1, "0"),), ".74", geometry,
                      ((1, 10, 16, 18), (6, 7, 10, 18)))
    return source(root, ((CENTER, "0"), (CENTER+1, ".0294")), ".1", (),
                  ((6, 7, 10, 13), (6, 7, 10, 12)))


def compare(bounds, windows):
    out = []
    for (lo, hi), (a, b) in zip(bounds, windows, strict=True):
        if hi < arb(a):
            out.append("低")
        elif lo > arb(b):
            out.append("高")
        else:
            out.append("未分离")
    return out


def verify(root):
    transport = run_path(str(root / "ch39-adjacent-pvm-transport.py"))
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    training_source = run_path(str(root / "ch39-common-embedding-central-right-ckm.py"))["verify_nodes"]
    four = run_path(str(root / "ch39-free-four-site-native-positive-spectrum.py"))
    reference = run_path(str(root / "ch39-three-shell-outer-holdout.py"))
    ck = run_path(str(root / "ch39-shell-clock-ckm-blocks.py"))
    flavor = run_path(str(root / "ch39-three-gram-source-flavor.py"))
    ctx.prec = 320
    angles = [arb(s).sin() for s in (".22704", ".04275", ".00377")]
    mixing = acb_mat(ck["standard_matrix"](*angles, arb("1.139")))
    ref_prob = [abs(z)**2 for z in mixing.entries()]
    ref_j = (mixing[0, 0]*mixing[1, 1]*mixing[0, 1].conjugate()*mixing[1, 0].conjugate()).imag
    print("参考 CKM 模平方及 J", [z.str(13) for z in (*ref_prob, ref_j)])
    for sign in (-1, 1):
        frame, trained, _ = training(root, sign, training_source)
        start = CENTER + sign
        bases = {}
        weights = {}
        for charge in CHARGES:
            _, bases[charge] = transport["spectrum"](shell, charge, start)
            weights[charge] = trained[start, charge]
        for distance in (2, 3):
            number = CENTER + sign*distance
            sectors = four["native"](root, number)
            projectors = []
            for i, (charge, native) in enumerate(zip(CHARGES, sectors, strict=True)):
                _, basis = transport["spectrum"](shell, charge, number)
                value = transport["transport"](bases[charge], basis, weights[charge])
                bases[charge], weights[charge] = basis, value
                assert all(0 < w < 1 for w in value)
                original = sum((w*p for w, p in zip(value, native, strict=True)), acb_mat(4, 4))
                bounds = reference["interlacing"](original)
                windows = reference["WINDOWS"][i]
                print("任意三平面", number, charge, [(lo.str(13), hi.str(13)) for lo, hi in bounds],
                      compare(bounds, windows), flush=True)
                effective = four["adj"](frame)*original*frame
                levels = sorted((v.real for v in effective.eig(algorithm="rump")), key=lambda v: v.mid())
                assert all(v > 0 for v in levels)
                print("同侧训练面", number, charge, [y.str(13) for y in levels],
                      compare([(y.lower(), y.upper()) for y in levels], windows), flush=True)
                if not all(x < y for x, y in zip(levels, levels[1:])):
                    continue
                squared, proj = flavor["projections"](effective*effective)
                assert all(z.overlaps(v*v) for z, v in zip(squared, levels, strict=True))
                projectors.append(proj)
            if len(projectors) != 2:
                print("固定面 CKM 未定", number)
                continue
            up, down = projectors
            prob = [(pu*pd).trace().real for pu in up for pd in down]
            j = (up[0]*down[0]*up[1]*down[1]).trace().imag
            ck_values = [*prob, j]
            fit = ["高" if y > t else "低" if y < t else "未分离"
                   for y, t in zip(ck_values, [*ref_prob, ref_j], strict=True)]
            print("固定面 CKM 模平方及带号 J", number, [z.str(13) for z in ck_values], fit, flush=True)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
