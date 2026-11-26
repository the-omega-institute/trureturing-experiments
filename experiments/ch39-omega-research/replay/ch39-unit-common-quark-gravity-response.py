"""一个全实轴介于零和一的有理响应同时容纳荷六截止与真实零点的夸克中心。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def top_feedback(rows, top_squared, cutoffs, pivot, heat, newton):
    # 破缺场表第二项是含三色和自旋计数的轻顶，z=m_t²/M_P²。
    multiplicity, coefficient, old_z = rows[1]
    scale = (arb("246.22")/arb("1.22089e19"))**2/2
    assert multiplicity == -12 and coefficient == 1
    assert old_z.overlaps(arb(".8616")**2*scale)
    new_z = top_squared*scale
    path = old_z.union(new_z)
    assert path > 0
    # d[x E₂(z/x)]/dz = −E₁(z/x)。用两质量之间的区间控制每层差，
    # 再以有符号首枢轴谱权求和，保留相消与所有参数误差。
    exposure = sum(c.real*heat(1,path/x) for c,x in zip(pivot,cutoffs,strict=True))
    delta = -(new_z-old_z)*coefficient*exposure/(2*arb.pi())
    corrected = newton+delta
    assert corrected < 0 and abs(delta) < arb("2.4e-33")
    return delta,corrected


def verify(root):
    shell = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    native = run_path(str(root / "ch39-free-four-site-native-positive-spectrum.py"))
    quarks = run_path(str(root / "ch39-complex-four-site-unit-quark-amplitudes.py"))
    gravity = run_path(str(root / "ch39-three-index-free-atoms.py"))
    pivots = run_path(str(root / "ch39-native-rh-ray-pvm.py"))
    ordered = run_path(str(root / "ch39-bounded-common-rational-amplitude.py"))["ordered"]
    energy = run_path(str(root / "ch39-three-gram-source-flavor.py"))["projections"]
    frame, weights = quarks["verify"](root)
    center = 2_306_888_100_085
    ctx.prec = 320
    height = acb.zeta_zero(center).imag
    data = []
    for charge,w in zip((-4,2),weights,strict=True):
        data.extend(zip(ordered(shell,charge,height),w,strict=True))
    # 取已由三条件共同确定的 Hill 根，以下各球是同一根的外包。
    h = arb("5.66600e-35").union(arb("5.66604e-35"))
    cap = arb(".635595244").union(arb(".635595246"))
    s = arb("3.291030650e-40").union(arb("3.291030654e-40"))
    six = ordered(shell,6,height)
    for i,x in enumerate(six):
        power = arb(1) if i==18 else (x/six[-1])**6000
        data.append((x,(h+cap*power/(power+s)).sqrt()))
    data.sort(key=lambda z:z[0].mid())
    assert len(data)==57 and all(v[0]<w[0] for v,w in zip(data,data[1:]))
    assert all(0<w<arb(".948") for t,w in data)

    slopes = [abs(w-v)/(t-x) for (x,v),(t,w) in zip(data,data[1:])]
    peak = max(range(len(slopes)), key=lambda i:slopes[i].mid())
    assert all(slopes[peak]>v for i,v in enumerate(slopes) if i!=peak)
    assert arb("116.0574429194") < slopes[peak] < arb("116.0574429196")
    print("最小 Lipschitz 常数、两个相邻中央谱点及其振幅",slopes[peak].str(13),
          *[v.str(12) for pair in data[peak:peak+2] for v in pair],flush=True)

    def response(x):
        factors = [1/(x-t)**2 for t,w in data]
        assert all(f>0 for f in factors)
        return sum(f*w for f,(t,w) in zip(factors,data,strict=True))/sum(factors)

    moments,heat,fields = gravity["broken_fields"](root,(center-1,center+1))
    for shift,rows in zip((-1,1),fields,strict=True):
        n = center+shift
        hh = acb.zeta_zero(n).imag
        projections = native["native"](root,n)
        for q,projectors in zip((-4,2),projections,strict=True):
            coeffs = [response(t) for t in ordered(shell,q,hh)]
            full = sum((w*p for w,p in zip(coeffs,projectors,strict=True)),acb_mat(4,4))
            m = native["adj"](frame)*full*frame
            values,_ = energy(m*m)
            ratio = (values[0]/values[-1]).sqrt()
            print("偏移、荷、最轻与最重夸克振幅比",shift,q,ratio.str(10),flush=True)
            if q==-4:
                top_squared = values[-1]
                assert ratio>arb(".00000629")/arb(".8573")
        _,levels,pivot,_ = next(t for t in pivots["shell_weights"](root,n) if t[0]==6)
        x = [response(t)**2 for t in levels]
        assert all(0<v<1 for v in x)
        columns = [moments(rows,heat,v) for v in x]
        K,W = [sum(c.real*pair[j] for c,pair in zip(pivot[0],columns,strict=True)) for j in range(2)]
        print("偏移、最小最大截止 Q/M_P、K、W",shift,min(v.lower() for v in x).sqrt().str(10),
              max(v.upper() for v in x).sqrt().str(10),K.str(12),W.str(12),flush=True)
        assert not K.contains(1) and not W.contains(arb("1e-122"))
        delta,corrected = top_feedback(rows,top_squared,x,pivot[0],heat,K)
        print("偏移、反馈顶 Yukawa、牛顿项改变量、反馈后 K",shift,
              top_squared.sqrt().str(12),delta.str(12),corrected.str(12),flush=True)


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
