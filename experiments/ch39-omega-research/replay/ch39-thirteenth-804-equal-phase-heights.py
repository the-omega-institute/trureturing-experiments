"""按登记次序在 804 三角相位相同的连续高度上重算完整支撑指标。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def verify(root):
    def lib(name):
        return run_path(str(root / ("ch39-" + name + ".py")))

    number = 2_306_888_100_098
    shell = lib("three-ray-chiral-boundary")
    geo = shell["geo"]
    gram = lib("native-rh-ray-pvm")["gram_cholesky"](root)[2]
    read = lib("eighth-ninth-wilson-height")["read_weights"]
    table = lib("cross-shell-pivot-support")["projected_tables"]
    core = lib("fourth-third-pivot-global-kinetic-obstruction")
    arc = lib("eighth-shortest-arc-pivot-support")
    ctx.prec = 320
    origin = acb.zeta_zero(number).imag
    period = 2*arb.pi()/(6*geo["b"])
    assert period > 0 and arc["SCALE"] == 10**7
    pole = shell["soft"]["pole"]
    i, j = next((i, j) for i, j, d, _ in geo["edges"] if d == 804 and pole not in (i, j))
    for k in range(9):
        height = origin + k*period
        phases = geo["phases"](height)
        magnetic = shell["laplacian"](6, height, phases=phases)
        loop = -magnetic[pole, i]*magnetic[i, j]*magnetic[j, pole]
        loop = loop/abs(loop)
        if k == 0:
            original_loop = loop
            print("三角顶点、起点 Wilson 相位和高度步长", (pole, i, j), loop.str(16),
                  period.str(18), flush=True)
        else:
            # 精确相等来自 6b(h_k-h_13)=2πk；这里直接核对原矩阵的有向闭路乘积。
            assert loop.overlaps(original_loop)
        weights = read(shell, gram, height, phases)
        tables = table(core, root, number, None, height=height)
        bounds = arc["score"](tables, weights)
        if bounds is None:
            assert k != 0
            print("k、高度、指标未定", k, height.str(29), flush=True)
            continue
        lo, hi = map(int, bounds)
        print("k、高度、S3 区间分子和分母", k, height.str(29), lo, hi, arc["SCALE"], flush=True)
        if k == 0:
            assert 0 < 13283 <= lo <= hi <= 13287
        elif hi < 0:
            return k, lo, hi
    return None


if __name__ == "__main__":
    result = verify(Path(__file__).resolve().parent)
    assert result is not None
    k, lo, hi = result
    assert k == 1 and -8687453 <= lo <= hi <= -8687449
