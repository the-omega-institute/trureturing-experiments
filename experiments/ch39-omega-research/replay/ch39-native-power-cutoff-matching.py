"""实际荷二壳与编号二的 Gram 枢轴的共同幂读出，同时标定牛顿项与弱差。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def native_readout(root):
    ctx.prec=256
    native=run_path(str(root/"ch39-native-rh-ray-pvm.py"))
    charge,levels,direct,mirror=next(item for item in native["shell_weights"](root) if item[0]==2)
    weights=[z.real for z in direct[1]]
    assert charge==2 and sum(weights).contains(1)
    assert all(w.overlaps(m.real) for w,m in zip(weights,mirror[1],strict=True))
    ratios=[e/levels[-1] for e in levels[:-1]]+[arb(1)]
    assert ratios[0]>0 and all(a<b for a,b in zip(ratios,ratios[1:]))
    _,fields=run_path(str(root/"ch39-three-index-free-atoms.py"))["spectra"](root)
    ctx.prec=256
    heat=run_path(str(root/"ch39-chiral-minimal-neutral-gravity.py"))["expint"]
    moment=run_path(str(root/"ch39-signed-cutoff-gravity.py"))["moments"]
    weak=native["weak_series"]()
    mp=arb("1.22089e19")
    p0=(arb("91.1876")/mp)**2
    baseline=arb("-.074256577194663365546").union(arb("-.074256577194663365544"))
    def moments(xs):
        values=[moment(fields[0],heat,x) for x in xs]
        return tuple(sum(c*v[i] for c,v in zip(weights,values,strict=True)) for i in range(2))

    def weak_difference(xs):
        x,*others=xs
        assert x>p0 and all(v>x for v in others)
        err=sum(abs(c)*(arb("2.398e-21")+arb("7.468e-36")/v)
                for c,v in zip(weights[1:],others,strict=True))
        return weights[0]*weak(x/p0)+(1-weights[0])*baseline+arb(0,err.upper())

    tv=sum(abs(c) for c in weights)
    massless=sum(a for _,a,_ in fields[0])/(2*arb.pi())
    # x(1-E2(z/x)) 在 0<x<=1 上递增。
    newton_error=tv*sum(abs(a)*(1-heat(2,z)) for _,a,z in fields[0])/(2*arb.pi())
    floor=(1-newton_error)/abs(massless*tv)
    low_degree_upper=baseline+tv*arb("2.398e-21")+arb("7.468e-36")*sum(
        abs(c)/ratio**12 for c,ratio in zip(weights,ratios,strict=True))/floor
    assert floor>arb(".009725") and newton_error<arb(".004901") and low_degree_upper<arb("-.06649")
    return dict(weights=weights, ratios=ratios, moments=moments, weak=weak_difference,
                mp=mp, p0=p0, baseline=baseline, massless=massless, newton_error=newton_error,
                floor=floor, low_degree_upper=low_degree_upper)


def verify(root):
    model=native_readout(root)
    weights,ratios,mp,p0,baseline=(model[k] for k in ("weights","ratios","mp","p0","baseline"))
    def cutoffs(scale,exponent):
        return [scale*ratio**exponent for ratio in ratios]
    def norm(scale,exponent):
        return model["moments"](cutoffs(scale,exponent))[0]
    def weak_difference(scale,exponent):
        return model["weak"](cutoffs(scale,exponent))

    left,right=arb(".02634069"),arb(".02634070")
    low,high=arb("12.883072"),arb("12.883074")
    scales,degrees=left.union(right),low.union(high)
    nleft,nright=norm(left,degrees)-1,norm(right,degrees)-1
    dlow,dhigh=weak_difference(scales,low),weak_difference(scales,high)
    print("整面 Newton 残差及整面弱差",*[v.str(13) for v in (nleft,nright,dlow,dhigh)],flush=True)
    assert nleft<0<nright and dlow<0<dhigh
    xs=cutoffs(scales,degrees)
    print("首末截止质量及第二层质量（GeV）",
          *[(mp*x.sqrt()).str(14) for x in (xs[0],xs[-1],xs[1])],flush=True)
    assert xs[0]>0 and scales<1
    vacuum=model["moments"](xs)[1]
    print("整盒无量纲真空读数 W",vacuum.lower().str(15),vacuum.upper().str(15),flush=True)
    assert arb("-.117248196") < vacuum < arb("-.117248059")
    floor,massless,newton_error,low_degree_upper=(model[k] for k in
        ("floor","massless","newton_error","low_degree_upper"))
    def allowed_scale(d):
        slope=massless*sum(c*s**d for c,s in zip(weights,ratios,strict=True))
        assert slope>0
        return (1+arb(0,newton_error.upper()))/slope
    middle=weak_difference(allowed_scale(arb(13)),arb(13))
    print("K 标定的全幂尺度下界、d<=12 的弱差上界、d=13 的弱差",
          floor.lower().str(12),low_degree_upper.upper().str(12),middle.str(12),flush=True)
    assert arb(".087")<middle<arb(".093")
    # 已证负的固有时密度给 D(x) 严格递减。所有非平凡尾权严格为负。
    tails=[-sum(weights[j:]) for j in range(1,len(weights))]
    sigma=min(v.lower() for v in tails)
    assert sigma>arb("3.1695")
    scale14=allowed_scale(arb(14).union(15))
    small=max((scale14.upper()*ratios[0]**14).upper(),(ratios[0]**15).upper())
    first_error=2*small*(1/p0+1/(arb("1.176397e14")/mp)**2)
    last_upper=baseline+arb("2.398e-21")+arb("7.468e-36")/floor
    positive=-sigma*last_upper-(1+sigma)*first_error
    print("d>=14 的最小尾负量、弱差下界",sigma.str(12),positive.lower().str(12),flush=True)
    assert positive>arb(".1718")


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
