"""同一荷在两个实际高度精确重建十九层的电导漂移下界。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def verify():
    root = Path(__file__).resolve().parent
    shared = run_path(str(root/"ch39-shared-even-heavy-towers.py"))
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    ctx.prec = 224
    n = 2_306_888_100_085
    data = {}
    for index in (11661,n,n+1):
        h = acb.zeta_zero(index).imag
        delta = arb(23)/25 * arb(2)**-20 / (h*h+shell["soft"]["point"]["study"]["offset"])
        c = 15-delta/3
        sums, fourth, moments = shared["target_data"](shell,h,c)
        assert all(t < 1500 for options in moments.values() for t,z in options)
        data[index] = sums, moments, delta

    for index, floors in ((11661,(".02034",".08022",".06567")),
                          (n+1,(".02589",".0012630",".0035214"))):
        cap = 1500/(1-data[n][2])*abs(data[n][2]-data[index][2])
        assert cap < arb("1e-8")
        for q, bound in zip((2,4,6),floors):
            a = data[n][1][q]
            b = data[index][1][q]
            gap = min(abs(v-w).lower() for v,_ in a for w,_ in b)
            assert gap > arb(bound)
            print("目标荷、比较编号、原始迹最小距离、总电导改动的必要下界",
                  q,index,gap.str(14),(gap/2).str(12))
        print("仅按原亏额更新时的原始迹变化上界",index,cap.str(12))
    # 放宽三个轻质量，只需它们各小于 10，留下至多 60 的两点比较余量。
    # 十六项精确重目标与三轻界仍给出 T<1500，三轻的最大原始根小于 16。
    sums0, moments, delta = data[n]
    assert all(t+3*16 < 1500 for t,z in moments[6])
    diff = abs(delta-data[n+1][2])
    total = 1500/(1-delta)
    # Q=wᵀHw，H 逐项非负；0≤w≤a 且两点 |Δw|≤a|Δδ|。
    # Q(a)≤(2∑a)²，对 cT 的变化另计 |Δc|T+c|ΔT|。
    drift = diff*(2*total**2+total/3+15*total)
    assert drift < arb("1e-15")
    variation = 32*abs(sums0[6]-data[n+1][0][6])
    assert variation-60-drift > arb("22.5054")
    print("荷六两处目标重层和差、扣除三轻及径向变化后的缺口",
          variation.str(16),(variation-60-drift).str(16))


if __name__ == "__main__":
    verify()
