"""逐编号排查第 12208 至 13700 枚实际零点，列出同时通过三代分界与两圈强耦合窗口的全部候选。

判据与 ch39-rotation-character-portal.py 的四枚规范候选相同：
spectral_supplier 的五种荷在分界 9 处各留三个低模，十六标量、钟弦分母 724 的两圈交点
反推 0.117 < alpha_s(M_Z) < 0.119。三模判据不成立时要求至少一种荷的第三、第四层被证书分在 9 的同侧；
三模成立时，先用 crossing_enclosure 同一像映射检查窗口种子盒 L∈29.040±0.014、1/alpha∈1/.1184±0.04：
像与盒不相交即认证“种子盒内无交点”，再用第二个种子 (29.166, .1212) 认证一个落在 (.1211, .1214) 的交点；
像与盒相交时 crossing_enclosure 必须成功并给出窗内强耦合，否则脚本报错。
结论只在“由种子包围的两圈交点分支”这一约定下成立（交点唯一性未证，约定同 ch39-chiral-clock-cohort.py）：
种子盒只覆盖 alpha∈(0.11784,0.11896)、L∈[28.986,29.094]，窗口其余部分与盒外分支没有排除。
判据 3（Q 以下恰有十四个荷三活跃层）不在本脚本检查，由 ch39-cubic-portal-holdout.py 中的 candidate 对五枚命中逐一断言；
“判据成立的恰为五枚”要靠两个脚本合起来成立。
"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, ctx

ROOT = Path(__file__).parent
fixed = run_path(str(ROOT / "ch39-chiral-clock-unit-readout.py"))
chiral = fixed["chiral"]

START, STOP = 12208, 13700
HITS = (12336, 12610, 12886, 13291, 13699)
OUTSIDE = (12409, 12546, 12675, 12813, 12950, 13088, 13089, 13227, 13495)


def seed_box_empty(model, initial, layers):
    """复写 crossing_enclosure 的像映射；像与窗口种子盒不相交即认证盒内无交点。"""
    two, sm = fixed["two"], fixed["sm"]
    unit = arb("2.435e18")/724**2
    start = (unit/model["mz"]).log()
    table = [(start+l, w, B, f, s) for l, w, B, f, s in two["thresholds"](model, layers, False)]
    L = arb("29.040")+arb(0, ".014")
    inverse = 1/arb(".1184")+arb(0, ".04")
    pi = arb.pi()
    definite = [r for r in table if r[0] < L]
    possible = [r for r in table if not (r[0] > L)]
    segments = fixed["gauge_segments"](L, possible, 8192)
    corr, _ = two["flow"]([*initial[:2], inverse], segments, sm["D"], sm["C"], min_inverse=8)
    rhs, mat = [], []
    for i in range(3):
        b0 = sm["BETA"][i]+sum(r[1][i] for r in definite)
        a = initial[i] if i != 2 else 0
        val = 2*pi*(a+corr[i])+sum(r[1][i]*r[0] for r in definite)
        val -= sum(w[i]*(L-l).max(0) for l, w, _, _, _ in possible[len(definite):])
        rhs.append([val])
        mat.append([b0, 2*pi, -2*pi*int(i == 2)])
    result = arb_mat(mat).solve(arb_mat(rhs)).entries()
    return not L.overlaps(result[0]) or not inverse.overlaps(result[2])


def classify(model, initial, read, height):
    rows = read(height)
    layers = {q: chiral["ordered"]([h for h, z in rows if z == q]) for q in (0, *model["charges"])}
    if not all(layers[q][2] < 9 < layers[q][3] for q in model["charges"]):
        assert any(layers[q][2] > 9 or layers[q][3] < 9 for q in model["charges"])
        return "no-three", None
    ctx.prec = 192
    if seed_box_empty(model, list(initial), layers):
        alpha = fixed["solve"](model, list(initial), layers, False, 724, "29.166", ".1212")
        assert arb(".1211") < alpha < arb(".1214")
        return "outside", alpha
    _, _, alpha = fixed["crossing_enclosure"](model, list(initial), layers, False, 724, "29.040", ".1184")
    assert arb(".117") < alpha < arb(".119")
    return "window", alpha


def verify():
    model, _, _, initial, _ = chiral["source"]["mean"]["inputs"]()
    ctx.prec = 128
    read = chiral["source"]["spectral_supplier"](model)
    hits, outside = [], []
    n = START
    while n <= STOP:
        count = min(250, STOP + 1 - n)
        ctx.prec = 128
        zeros = acb.zeta_zeros(n, count)
        assert len(zeros) == count
        for k, zero in enumerate(zeros):
            ctx.prec = 128
            kind, alpha = classify(model, initial, read, zero.imag)
            if kind == "window":
                hits.append(n + k)
                print("规范窗口候选", n + k, alpha.str(8), flush=True)
            elif kind == "outside":
                outside.append(n + k)
                print("三模、窗口种子盒内无交点、第二种子交点强耦合", n + k, alpha.str(6), flush=True)
        n += count
    assert tuple(hits) == HITS, hits
    assert tuple(outside) == OUTSIDE, outside
    print("第", START, "至", STOP, "枚中规范窗口候选", hits, "三模窗外", outside)


if __name__ == "__main__":
    verify()
