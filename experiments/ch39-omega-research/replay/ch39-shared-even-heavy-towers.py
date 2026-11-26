"""三种偶荷共用电导时，二步质量的有符号迹限制跨荷重塔。"""
from itertools import combinations, product
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def target_data(shell, h, c):
    sums, fourth, moments = {}, {}, {}
    for q in (2,4,6):
        eig = shell["laplacian"](q,h).eig(multiple=True)
        assert len(eig) == 19 and all(z.imag.contains(0) and z.real > 0 for z in eig)
        m = sorted((z.real*abs(z.real-c) for z in eig),key=float)
        assert all(x < y for x,y in zip(m,m[1:]))
        assert m[2] < 10 < m[3]
        sums[q], fourth[q] = sum(m[3:]), m[3]
        upper = sum((c+(c*c+128*z).sqrt())/2 for z in m[3:])
        assert 32*m[3] > c*c/4 and 4*m[2] < c*c
        choices = [[(c+(c*c+4*z).sqrt())/2,(c+(c*c-4*z).sqrt())/2,
                    (c-(c*c-4*z).sqrt())/2] for z in m[:3]]
        heavy_square = sum(((c+(c*c+128*z).sqrt())/2)**2 for z in m[3:])
        moments[q] = [(upper+sum(choice), heavy_square+sum(z*z for z in choice))
                      for choice in product(*choices)]
    return sums, fourth, moments


def verify():
    root = Path(__file__).resolve().parent
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    ctx.prec = 224
    h = acb.zeta_zero(2_306_888_100_085).imag
    delta = arb(23)/25*arb(2)**-20 / (h*h+shell["soft"]["point"]["study"]["offset"])
    c = 15-delta/3
    sums, fourth, moments = target_data(shell,h,c)
    traces = {q: [t for t,z in options] for q, options in moments.items()}
    neutral = shell["laplacian"](0,h)
    t0 = sum(neutral[j,j].real for j in range(19))
    square = neutral**2
    q0 = sum(square[j,j].real for j in range(19))
    # 若所有电导只沿固定中性电导的一条直线变化，第二矩对第一矩至多为二次。
    slopes = {}
    for charge, options in moments.items():
        assert all(t > t0 for t,z in options)
        slopes[charge] = [(t-t0,(z-q0)/(t-t0)) for t,z in options]
    det_floor = min(abs((b-a)*(vc-va)-(e-a)*(vb-va)).lower()
                    for (a,va),(b,vb),(e,vc) in product(slopes[2],slopes[4],slopes[6]))
    print("原中性前两矩、单一电导变化方向的直线行列式绝对值下界",
          t0.str(14),q0.str(14),det_floor.str(15))
    assert det_floor > arb(".00002794")
    # 三个偶荷的谱迹必须满足逐边处方在 q² 上产生的同一插值关系。
    single = min(abs(5*a-8*b+3*z).lower()
                 for a,b,z in product(traces[2],traces[4],traces[6]))
    def divided_trace(trace, charge):
        t = charge**2
        return (trace-18*c)/(t*(t-1)*(t-9))
    fixed = [[divided_trace(v,q) for v in traces[q]] for q in (2,4,6)]
    anchored = min(abs(5*a-8*b+3*z).lower() for a,b,z in product(*fixed))
    assert single > arb(".00085735") and anchored > arb("43.25")
    print("三区谱迹的一次插值残差、固定中性及两奇荷后的四次插值残差",
          single.str(16),anchored.str(16))
    for q,r in combinations(sums,2):
        deficit = 32*abs(sums[q]-sums[r])
        assert deficit > 60
        trace_gap = min(abs(v-w).lower() for v in traces[q] for w in traces[r])
        assert trace_gap > arb(".39")
        print("两个偶荷、目标重层和差、保三轻质量的最小原始迹差",
              q,r,deficit.str(14),trace_gap.str(14))
    # 用整数剩余类直接给出荷依赖电导表的最小周期。
    def residue_class(q, period):
        r = q % period
        return min(r,period-r)
    for period in range(2,13,2):
        classes = {residue_class(q,period) for q in (0,2,4,6)}
        assert (len(classes) == 4) == (period == 12)
    print("电导依赖整数荷、保共轭且奇荷原权时的最小周期",12)
    gap = sums[6]-sums[2]
    kcap = 60/gap
    assert min(fourth.values())*kcap > c*c/4
    absolute = (32*gap-60)/32
    relative = (32*gap-60)/(32*(sums[6]+sums[2]))
    # 若尝试达到此误差以下，所有候选重层都处于正支。
    assert 32*min(fourth.values())-absolute > c*c/4
    assert 32*min(fourth.values())*(1-relative) > c*c/4
    assert kcap < arb("4.403") and absolute > arb("11.75") and relative > arb(".00730")
    print("共用电导的倍率必要上界、32倍逐层最大绝对与相对误差下界",
          kcap.str(15),absolute.str(15),relative.str(15))


if __name__ == "__main__":
    verify()
