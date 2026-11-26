"""由荷二重模的曲率窗在全部两钟角上排除中间相位，并约束负相位的味份额。"""
from functools import cache
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

root = Path(__file__).parent
window = run_path(str(root / "ch39-common-step-zero-window.py"))
cp = run_path(str(root / "ch39-cp-ensemble-test.py"))


def exclude_middle():
    pi = arb.pi()
    left = arb(7)/64
    span = 1-2*left

    @cache
    def count_at(depth):
        # 单元半宽为 π/(2·2^depth) 和 span 倍该数。
        norm = 4*(1+span)*pi/(2**depth)
        return window["supplier"](2,arb("18.94")-norm)

    @cache
    def center(depth, index, relative):
        angle = arb(2*index+1)/2**(depth+1)
        if relative:
            angle = left+span*angle
        return acb(0,pi*angle).exp()

    pending = [(3,i,j) for i in range(8) for j in range(8)]
    leaves = []
    while pending:
        depth, i, j = pending.pop()
        w = center(depth,j,False)
        z = w*center(depth,i,True)
        if count_at(depth)({724:z,804:w}) >= 15:
            leaves.append(depth)
            continue
        assert depth < 12
        pending.extend((depth+1,2*i+a,2*j+b) for a in range(2) for b in range(2))
    assert len(leaves) == 10744 and max(leaves) == 11
    assert sum(4**(11-k) for k in leaves) == 4**11
    print("覆盖整个中间相位带的已证单元数与最大二分深度", len(leaves), max(leaves))


def verify():
    ctx.prec, ctx.cap = 192, 3
    # 若荷二的第十六质量落在三模、曲率两窗，先用大于 9 的分界控制截止余项。
    mass_unit = 32*arb("2.435e18")/724**2
    vev = arb("246.22")
    planck = arb("1.22089e19")
    coeff = (9*arb(3).log()-12*arb(2).log())/2
    low_sector_max = (vev**2/(9*mass_unit*planck))**2
    fraction = 1-low_sector_max/(arb(".03")*coeff)
    mass_floor = vev**2/(mass_unit*planck)*(fraction*coeff/(arb.pi()*arb("1.002e-122")))**(arb(1)/4)
    assert mass_floor > arb("74.64")
    # 中性质量中心是 15−δ/3，软化半径权只减小磁拉普拉斯。
    perturbation = arb(23)/25*arb(2)**-20/14**2
    cutoff = arb("18.94")
    assert max(arb(225)/4, cutoff*(cutoff-15+perturbation/3)) < mass_floor
    exclude_middle()

    ctx.prec, ctx.cap = 384, 3
    flavor = cp["flavor"]
    axes = [flavor["singlet_at"]({724:z,804:acb(1)})[1] for z in (acb(1), acb(0,1))]
    lower = axes[0]*(7*arb.pi()/64).cos()-abs(axes[1])
    upper = sum(t*t for t in axes).sqrt()
    assert lower > arb(".9404240")
    data = cp["two"]["spectral_data"]()
    for j in range(128):
        c = (arb(-1)+arb(j)/64).union(-1+arb(j+1)/64)
        assert cp["two"]["probabilities"](c,data)[0][1] < 0
    def probability(c):
        return cp["two"]["probabilities"](c,data)[0][0]
    up0 = probability(flavor["singlet"](11661)[1])
    min_up = probability(upper)
    neg_floor = probability(-lower)
    fraction = (up0-min_up)/(neg_floor-min_up)
    assert neg_floor > arb("127.32")*up0 and fraction < arb(".00575853")
    print("味坐标绝对值下界、负支跃迁概率倍率下界与味份额上界",
          lower.str(13), (neg_floor/up0).str(13), fraction.str(13))


if __name__ == "__main__":
    verify()
