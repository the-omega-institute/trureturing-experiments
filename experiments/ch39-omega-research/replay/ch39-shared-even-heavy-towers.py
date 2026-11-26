"""三种偶荷共用电导时，二步质量的有符号迹限制跨荷重塔。"""
from itertools import combinations, product
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def verify():
    root = Path(__file__).resolve().parent
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    ctx.prec = 224
    h = acb.zeta_zero(2_306_888_100_085).imag
    delta = arb(23)/25*arb(2)**-20 / (h*h+shell["soft"]["point"]["study"]["offset"])
    c = 15-delta/3
    sums, traces, fourth = {}, {}, {}
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
        traces[q] = [upper+sum(choice) for choice in product(*choices)]
    for q,r in combinations(sums,2):
        deficit = 32*abs(sums[q]-sums[r])
        assert deficit > 60
        trace_gap = min(abs(v-w).lower() for v in traces[q] for w in traces[r])
        assert trace_gap > arb(".39")
        print("两个偶荷、目标重层和差、保三轻质量的最小原始迹差",
              q,r,deficit.str(14),trace_gap.str(14))
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
