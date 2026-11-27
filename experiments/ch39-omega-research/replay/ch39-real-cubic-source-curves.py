"""用两个有理区间界定位实系数三次厄米源的完整相位解集。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat


def verify():
    orbit = run_path(str(Path(__file__).with_name("ch39-two-zero-orbit-condition.py")))
    quadratic = run_path(str(Path(__file__).with_name("ch39-two-zero-quadratic-source.py")))
    model, up_map, down_map = orbit["calibration"]()
    adj, trace = model["adj"], quadratic["trace"]
    up, down = [adj(matrix) * matrix for matrix in (up_map, down_map)]
    up, down = up / trace(up).real, down / trace(down).real
    frame, _ = quadratic["normal_frame"](up, down)
    original = model["K"] * adj(model["K"])
    original_pairing = trace(frame[0] * original).real
    assert arb("-.0254750") < original_pairing < arb("-.0254748")
    cubics = (up * down * up, down * up * down)
    transverse = arb_mat([[trace(axis * word).real for word in cubics] for axis in frame[1:]])
    determinant = transverse.det()
    assert determinant < arb("-1.29e-8")

    commutator = up * down - down * up
    quartic = up * down * down * up
    pairing = trace(acb(0, 1) * commutator * quartic)
    determinant_pairing = -acb(0, 1) * commutator.det()
    assert pairing.real < arb("-2.32e-8") and determinant_pairing.real < arb("-2.32e-8")

    parts, _ = orbit["orbit"](model, up_map, down_map)
    harmonics = quadratic["orbit_coefficients"](parts, adj, frame[:1])
    constant, cross_cos, cross_sin, x_cos, x_sin, y_cos, y_sin = harmonics.entries()
    baseline = abs(constant) + (x_cos**2 + x_sin**2).sqrt()
    cross = (cross_cos**2 + cross_sin**2).sqrt()
    amplitude = (y_cos**2 + y_sin**2).sqrt() - cross
    assert baseline < arb(".029") and amplitude > arb(".055")
    # 在两个固定纬度，正弦项保证相反的严格符号。
    sign_margin = y_sin - cross - baseline
    assert sign_margin > arb(".027")
    slope = (amplitude**2 - baseline**2).sqrt()
    assert slope > arb(".048")
    print("三次增量的二阶子式", determinant.str(10))
    print("原协方差对三次超平面单位法向的迹配对", original_pairing.str(10))
    print("四次词对交换子方向的迹配对", pairing.real.str(10))
    print("常数项的模上界、关于 y 的振幅下界", baseline.upper().str(10), amplitude.lower().str(10))
    print("两固定纬度符号裕量与每个根的纵向斜率下界", sign_margin.lower().str(10), slope.lower().str(10))


if __name__ == "__main__":
    verify()
