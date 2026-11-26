"""固定荷六首枢轴 Hill 谱函数，在两枚相邻实际零点重算十九项。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def verify(root):
    center = 2_306_888_100_085
    numbers = (center,center-1,center+1)
    native = run_path(str(root/"ch39-native-rh-ray-pvm.py"))
    shell_data = []
    for number in numbers:
        _, levels, pivots, _ = next(entry for entry in native["shell_weights"](root,number)
                                     if entry[0] == 6)
        weights = [z.real for z in pivots[0]]
        assert len(levels) == len(weights) == 19 and sum(weights).contains(1)
        shell_data.append((levels,weights))
    moments,heat,fields = run_path(str(root/"ch39-three-index-free-atoms.py"))["broken_fields"](root,numbers)
    weak = run_path(str(root/"ch39-full-neutral-free-cutoff.py"))["weak_column"]
    ctx.prec = 256
    weak_low = native["weak_series"]()
    h = arb("5.66600e-35").union(arb("5.66604e-35"))
    kappa = arb(".635595244").union(arb(".635595246"))
    s = arb("3.291030650e-40").union(arb("3.291030654e-40"))
    mp = arb("1.22089e19")
    p0,p1 = [(arb(m)/mp)**2 for m in ("91.1876","1.176397e14")]
    target = (1,arb("1e-122"),0)
    holdouts = []
    eighteenth = ((".221603446",".221603448"),(".011962311",".011962313"),(".63250287",".63250289"))
    for number, (levels,weights), rows, bounds in zip(numbers,shell_data,fields,eighteenth,strict=True):
        powers = [(v/levels[-1])**6000 for v in levels[:-1]]+[arb(1)]
        sites = [h+kappa*r/(r+s) for r in powers]
        assert len(sites) == 19 and all(p0 < x < 1 for x in sites)
        assert arb(bounds[0]) < sites[17] < arb(bounds[1])
        columns = []
        for x in sites:
            if x < arb("1e-24"):
                # 已证 D(x) 递减，可用两端包住整个低截止区间。
                d = weak_low(x.lower()/p0).union(weak_low(x.upper()/p0))
            else:
                d = weak(x,p0,p1)
            columns.append((*moments(rows,heat,x),d))
        values = [sum(c*col[i] for c,col in zip(weights,columns,strict=True)) for i in range(3)]
        print("编号、第 16 至 19 层平方位置、K W D",number,
              [x.str(13) for x in sites[15:]],[v.str(15) for v in values],flush=True)
        residual = [v-g for v,g in zip(values,target,strict=True)]
        if number == center:
            assert all(r.contains(0) for r in residual)
        else:
            holdouts.append(residual)
    (km,wm,dm),(kp,wp,dp) = holdouts
    assert arb("-.46623104") < 1+km < arb("-.46623102")
    assert arb("1.53675837") < target[1]+wm < arb("1.53675840")
    assert arb("-.020727") < dm < arb("-.020725")
    assert arb("4.53372107") < 1+kp < arb("4.53372111")
    assert arb("-12.9389126") < target[1]+wp < arb("-12.9389124")
    assert arb(".010167") < dp < arb(".010169")
    # 相邻点总质量取 1。两项引力残差单调方向相反，交点处都为负。
    assert arb(".29324") < -km/(kp-km) < arb(".29326")
    assert arb(".10615") < -wm/(wp-wm) < arb(".10617")
    fraction = (wm-km)/(kp-km-wp+wm)
    newton = km+fraction*(kp-km)
    vacuum = wm+fraction*(wp-wm)
    weak_mix = dm+fraction*(dp-dm)
    assert 0 < fraction < 1 and newton < 0 and newton.overlaps(vacuum)
    assert abs(weak_mix) < -newton
    assert arb(".154192210") < fraction < arb(".154192216")
    assert arb(".69527732") < -newton < arb(".69527738")
    print("右邻最优份额、三个残差的最小最大绝对值和该处弱差",
          fraction.str(14),(-newton).str(14),weak_mix.str(13))



if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
