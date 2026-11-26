"""在实际上下相邻零点重建场及 PVM，检验中心荷六首枢轴三段解。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def group_weights(shell_weights, root, number):
    _, levels, direct, _ = next(entry for entry in shell_weights(root,number) if entry[0] == 6)
    assert len(levels) == 19 and all(x<y for x,y in zip(levels,levels[1:]))
    weights = [c.real for c in direct[0]]
    grouped = [sum(weights[:17]), weights[17], weights[18]]
    assert sum(grouped).contains(1)
    return grouped


def verify(root):
    n = 2_306_888_100_085
    numbers = (n,n-1,n+1)
    native = run_path(str(root/"ch39-native-rh-ray-pvm.py"))
    groups = [group_weights(native["shell_weights"],root,number) for number in numbers]
    gravity = run_path(str(root/"ch39-three-index-free-atoms.py"))
    moments,heat,fields = gravity["broken_fields"](root,numbers)
    ctx.prec = 256
    weak = run_path(str(root/"ch39-full-neutral-free-cutoff.py"))["weak_column"]
    weak_low = native["weak_series"]()
    high_pair = run_path(str(root/"ch39-ordered-cross-gravity-roots.py"))["high_pair"]
    mz, mp = arb("91.1876"), arb("1.22089e19")
    p0, p1 = [(m/mp)**2 for m in (mz,arb("1.176397e14"))]
    q_lo, q_hi = arb("97.16194"), arb("97.16196")
    t_lo, t_hi = (q_lo/mz)**2, (q_hi/mz)**2
    low = p0*t_lo.union(t_hi)
    middle, high = high_pair(groups[0], low, fields[0], heat, moments)
    a0, b0, c0 = groups[0]
    dm, dh = weak(middle,p0,p1), weak(high,p0,p1)
    for t, sign in ((t_lo,1),(t_hi,-1)):
        assert sign*(a0*weak_low(t)+b0*dm+c0*dh) > 0
    assert arb(".2215232227") < middle < arb(".2215232229")
    assert arb(".6353652947") < high < arb(".6353652949")
    print("中心三段物理截止，最低 GeV、其余两个平方位置",
          q_lo.str(9),q_hi.str(9),middle.str(14),high.str(14))
    # 在中心精确三式解处，低层弱差由另两层给出。三编号的轻场参数相同。
    dl = -(b0*dm+c0*dh)/a0
    assert dl.overlaps(weak_low(t_lo.union(t_hi)))
    target = [1,arb("1e-122"),0]
    residuals = []
    for number, weights, rows in zip(numbers,groups,fields,strict=True):
        cols = [(*moments(rows,heat,x),d) for x,d in zip((low,middle,high),(dl,dm,dh),strict=True)]
        residual = [sum(c*col[i] for c,col in zip(weights,cols,strict=True))-y
                    for i,y in enumerate(target)]
        print("编号、三段总权、K W D 残差",number,[v.str(11) for v in weights],
              [v.str(15) for v in residual],flush=True)
        if number == n:
            assert all(val.contains(0) for val in residual)
        else:
            residuals.append(residual)

    (km,wm,dmix),(kp,wp,dpix) = residuals
    assert arb("-.226981966") < km < arb("-.226981963")
    assert arb(".229547114") < wm < arb(".229547116")
    assert arb("-.008592883219") < dmix < arb("-.008592883217")
    assert arb(".242339265") < kp < arb(".242339267")
    assert arb("-.249123374") < wp < arb("-.249123372")
    assert arb(".009120657104") < dpix < arb(".009120657107")
    # 右邻在邻点总质量中的比例 t。K(t) 递增，W(t) 递减。
    t = (wm-km)/(kp-km-wp+wm)
    k = km+t*(kp-km)
    w = wm+t*(wp-wm)
    d = dmix+t*(dpix-dmix)
    assert 0 < t < 1 and k < 0 and k.overlaps(w) and abs(d) < -k
    print("最优邻点混合的右邻份额、三个残差的最大绝对值、该处弱差",
          t.str(14),(-k).str(14),d.str(14))
    assert arb(".481574965") < t < arb(".481574967")
    assert arb(".000968608") < -k < arb(".000968610")


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
