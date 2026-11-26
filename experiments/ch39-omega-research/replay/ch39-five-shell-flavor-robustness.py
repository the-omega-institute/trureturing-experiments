"""六个 Yukawa 与四个 CKM 参数同时变化时，核验五壳的三空间下界。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, ctx


FRACTION = "1e-6"


def input_model(root, fraction):
    model = run_path(str(root / "ch39-global-independent-weak-frame-separation.py"))
    data = run_path(str(root / "ch39-four-site-physical-cap-window.py"))
    conventions = run_path(str(root / "ch39-shell-clock-ckm-blocks.py"))
    ctx.prec = 256
    assert fraction > 0
    inputs = [arb(c) + arb(0, (fraction * arb(s)).upper())
              for c, s in zip(data["CENTERS"], data["WIDTHS"], strict=True)]
    masses = (inputs[:3], inputs[3:6])
    assert all(0 < u and u < v and v < w and w < 1 for u, v, w in masses)
    assert all(0 < angle and angle < arb.pi() / 2 for angle in inputs[6:9])
    mixing = acb_mat(conventions["standard_matrix"](
        *(angle.sin() for angle in inputs[6:9]), inputs[9]))
    assert all(abs(z) > 0 for z in mixing.entries())
    quartet = mixing[0, 0] * mixing[1, 1] * (mixing[0, 1] * mixing[1, 0]).conjugate()
    assert quartet.imag > 0
    geom = model["geom"]
    up, down = [geom["diag"]([m / row[2] for m in row[:2]] + [arb(1)]) for row in masses]
    A, B = geom["adj"](mixing) * up * mixing, down
    assert not (A * B - B * A).det().contains(0)

    def polygon(blocks, chart, t):
        # 原裁剪函数仅通过 -t m_local,3/m_remote,3 使用铅笔参数。
        remote = 1 - chart
        nominal = model["mass"]
        adjusted = t * (masses[chart][2] / masses[remote][2]) / (nominal[chart][2] / nominal[remote][2])
        return model["polygon"](blocks, chart, adjusted)

    return dict(model, mass=masses, V=mixing, A=A, B=B, polygon=polygon)


def verify(root, fraction=FRACTION):
    five = run_path(str(root / "ch39-five-shell-oriented-axis-triangle.py"))
    outer = run_path(str(root / "ch39-left-right-independent-weak-frame-separation.py"))
    model = input_model(root, arb(fraction))
    model["azimuth"]["compile_cover"](root)
    rays, rows, frames = {}, {}, {}
    for offset in five["LAYERS"]:
        rays[offset], rows[offset], frames[offset] = five["prepare"](model, offset)
    for chart in range(2):
        for branch in range(3):
            five["cover_branch"](model, outer, chart, branch, rows, frames)

    for i, j in ((-2, -1), (-1, 0), (-2, 0)):
        special = {(7, 7), (10, 7)} if (i, j) == (-2, 0) else set()
        distances = {(k, l): sum((a - b) ** 2 for a, b in zip(x, y, strict=True)).sqrt()
                     for k, _, x in rays[i] if k in five["LAYERS"][i]
                     for l, _, y in rays[j] if l in five["LAYERS"][j]}
        assert min(d.lower() for pair, d in distances.items() if pair not in special) > 2 * five["CHORD"]
        assert all(distances[pair] < 2 * five["CHORD"] for pair in special)
    assert arb("-.30") > arb("-.32")
    print("五壳三空间下界的十维联合盒半宽倍数", fraction)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
