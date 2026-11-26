"""相邻两个实际高度共用正截止、非负轮廓时，同时精确匹配引力两项。"""
from math import factorial
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, ctx


def verify():
    root = Path(__file__).resolve().parent
    previous = run_path(str(root / "ch39-two-zero-cutoff-gravity.py"))
    single, lib, rows0, rows1 = previous["adjacent_fields"](root)
    ctx.prec = 768
    epsilon = arb("1e-122")
    cutoff = epsilon/93
    width = arb("1e-20")
    small = cutoff*arb("1e-35")
    t1, t2 = cutoff*(1-width), cutoff*(1+width)
    peak_unit = 1116*arb.pi()/epsilon
    high_unit = arb("1e-200")
    def read(rows, x):
        if x > 2*cutoff:
            return single["moments"](rows, lib["expint"], x)
        # E_j(v)≤e^{-v}/(j-1)，用幂界控制红外重质量尾；轻层保留精确的二次式。
        k, w = x/(12*arb.pi()), 31*x*x/(4*arb.pi())
        zmin = min(z.lower() for _, _, z in rows[1:])
        assert zmin > arb("1e-11")
        bound = factorial(6)*(x/zmin)**6
        dk = x*bound*sum(abs(a) for _, a, _ in rows[1:])/(2*arb.pi())
        dw = x*x*bound*sum(abs(g) for g, _, _ in rows[1:])/(8*arb.pi())
        return k+arb(0, dk), w+arb(0, dw)
    def slab(left, right, scale):
        # 差矩另作归一。红外两背景公共的无质量项应在区间运算前消去。
        values = []
        for rows in (rows0, rows1):
            a = read(rows, left)
            b = read(rows, right)
            values.append(tuple((y-x)*scale for x, y in zip(a, b)))
        (k, w), (k1, w1) = values
        if right < 2*cutoff:
            bound = (factorial(6)*(right/arb("1e-11"))**6)*right*scale
            # 两背景各重层的绝对曲率权与偏振权之和均小于2000。
            assert all(sum(abs(a)+abs(g) for g,a,z in rows[1:]) < 2000 for rows in (rows0, rows1))
            dk, dw = arb(0, 4000*bound), arb(0, 4000*bound*right)
        else:
            dk, dw = k1-k, w1-w
        return [k, w/epsilon, dk/high_unit, dw/high_unit]
    base = slab(arb(0), small, arb(1))
    low = [slab(small,t1,peak_unit), slab(t1,t2,peak_unit)]
    mids = [arb("2.948e-11"), arb("1e-10"), arb("8.39532e-9")]
    relative_width = arb("1e-7")
    assert mids[-1] * (1 + relative_width) < arb("8.39533e-9")
    high = [slab(x, x*(1 + relative_width), high_unit) for x in mids]
    # 第一个中等截止窗的非负高度固定，剩下四个高度由两背景的四个目标决定。
    matrix = arb_mat([*low, *high[1:]]).transpose()
    assert not matrix.det().contains(0)
    rhs = arb_mat([[v-a-b] for v, a, b in zip((1,1,0,0), base, high[0])])
    sol = matrix.solve(rhs).entries()
    print("两个红外台阶的归一高度、后两个中间窗归一高度",*[v.str(18) for v in sol])
    assert all(v > 0 for v in sol)
    assert sol[0] < 1+arb("2e-20") and sol[1]<sol[0]
    assert all(v*high_unit < 1 for v in sol[2:])


if __name__ == "__main__":
    verify()
