"""一个全实轴介于零和一的有理响应同时容纳荷六截止与真实零点的夸克中心。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


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
                assert ratio>arb(".00000629")/arb(".8573")
        _,levels,pivot,_ = next(t for t in pivots["shell_weights"](root,n) if t[0]==6)
        x = [response(t)**2 for t in levels]
        assert all(0<v<1 for v in x)
        columns = [moments(rows,heat,v) for v in x]
        K,W = [sum(c.real*pair[j] for c,pair in zip(pivot[0],columns,strict=True)) for j in range(2)]
        print("偏移、最小最大截止 Q/M_P、K、W",shift,min(v.lower() for v in x).sqrt().str(10),
              max(v.upper() for v in x).sqrt().str(10),K.str(12),W.str(12),flush=True)
        assert not K.contains(1) and not W.contains(arb("1e-122"))


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
