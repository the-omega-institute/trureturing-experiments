"""直接用五对高斯有理圆点指定读弱的荷，验证整个弦长窗口。"""
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


CENTERS = ((1,2,21864900,-23471596), (-4,1,-35552782,46976248),
           (2,2,28741423,-33197008), (-3,1,33066878,-40836094),
           (6,2,-38552575,56752117))


def cayley(t):
    return acb(1-t*t,2*t)/(1+t*t)


def verify(root):
    graph = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    pump = run_path(str(root/"ch39-native-light-projector-pump.py"))
    ctx.prec = 160
    adjoint = pump["adjoint"]
    frames = list(pump["rotation_bases"](graph))
    heights = dict(run_path(str(root/"ch39-first-ray-forecast-phases.py"))["FORECAST"])
    tolerance = arb(1)/500
    limit = arb(".001000001")
    assert limit**2*(4-tolerance**2) > tolerance**2
    divisions = 16
    points = [cayley(limit*(2*j+1-divisions)/divisions) for j in range(divisions)]
    radius = 2*limit/divisions
    change = 2*arb(2).sqrt()*radius + 6*arb("1e-8")
    for q,k,a,b in CENTERS:
        frame = frames[k]
        index, = [j for j in range(frame.ncols()) if abs(frame[6,j]) > arb(".5")]
        frozen,links = pump["clock_matrices"](graph,q)
        base = adjoint(frame)*frozen*frame
        clocks = []
        h = acb.zeta_zero(heights[q]).imag
        theta = graph["geo"]["theta"]
        for link,integer,distance in zip(links,(a,b),(724,804),strict=True):
            anchor = cayley(arb(integer)/10**8)
            observed = acb(0,q*(theta(distance)-theta(244))*h).exp()
            assert abs(observed-anchor) < tolerance
            block = adjoint(frame)*link*frame
            terms = [block*(anchor*z) for z in points]
            clocks.append([term+adjoint(term) for term in terms])
        worst = arb(0)
        for first,second in product(*clocks):
            vector,gap = pump["first_line"](base+first+second)
            assert gap > 2*change
            sine = change/(gap-change)
            # 轨道坐标的模至多增加最低线主角的正弦。
            bound = abs(vector[index])+sine
            worst = max(worst,bound.upper())
        assert worst < arb(".000487")
        print("唯一读弱荷、指定字符、两圆盘内奇异值统一上界",q,k,worst.str(12))
    assert (tolerance/arb.pi())**2 > arb("4.05e-7")


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
