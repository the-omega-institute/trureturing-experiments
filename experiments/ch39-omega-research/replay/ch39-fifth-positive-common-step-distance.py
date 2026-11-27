"""第五壳五式与全程正动能使第四步远离原完整味见证。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb


def verify(root):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    center, radius, _, source, levels = lib("fifth-third-pivot-joint")["verify"](
        root, return_transport=True
    )
    shell = lib("three-ray-chiral-boundary")
    _, spectrum = lib("adjacent-pvm-transport")["spectrum"](shell, 6, 2_306_888_100_090)
    low, h, k = levels
    # 两列组合取同一组解析系数，允许各联合根的原上游精确值在源盒中变化。
    angle, phase = arb(".04357"), arb(".01451")
    v, w = angle.cos(), angle.sin() * acb(0, phase).exp()
    outgoing = spectrum.conjugate().transpose() * center**3 * source * acb_mat([[v], [w]])
    error = 3 * radius
    cap = (arb(".001") - low) * (abs(v)**2 / (h-low) + abs(w)**2 / (k-low))
    # 分别观测第 16、18、19 层，第四个结果是它们的整个正交补。
    prob = []
    for i in (15, 17, 18):
        mag = abs(outgoing[i, 0])
        assert mag.lower() > error
        prob.append((mag + arb(0, error.upper()))**2)
    remaining = 1 - sum(prob)
    lam = ((1 - 3*cap) / remaining).sqrt()
    assert min(lam**2 * u for u in prob) > cap
    # 在三个单层概率各自不超过 cap 的单纯形，亲和度最大值在三面交点。
    affinity = sum((u * cap).sqrt() for u in prob) + (remaining*(1-3*cap)).sqrt()
    separation = 2 * (affinity.acos()/6).sin()
    print("三项旧概率与新单项上限", *[u.str(15) for u in prob], cap.str(15), flush=True)
    print("最大亲和度、第四步射影算子距离下界", affinity.str(15), separation.str(15), flush=True)
    assert separation > arb(".29436")
    return separation


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
