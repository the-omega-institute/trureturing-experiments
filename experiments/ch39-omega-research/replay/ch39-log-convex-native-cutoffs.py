"""双对数凸的正保序谱标尺对中心壳三个联合条件的限制。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


def verify(root):
    native=run_path(str(root/"ch39-native-rh-ray-pvm.py"))
    fields=run_path(str(root/"ch39-three-index-free-atoms.py"))["broken_fields"]
    _,heat,(rows,)=fields(root,[2_306_888_100_085])
    ctx.prec=256
    alpha=-sum(a for _,a,_ in rows)/(2*arb.pi())
    error=sum(abs(a)*(1-heat(2,z)) for _,a,z in rows)/(2*arb.pi())
    mp=arb("1.22089e19")
    bandwidth=arb("2662.44")/arb("3.4741e16")
    p0=(arb("91.1876")/mp)**2
    d0=native["weak_series"]()(arb(1))
    base=arb("-.074256577194663365546").union(arb("-.074256577194663365544"))
    for q,levels,source,_ in native["shell_weights"](root):
        if q not in (1,-3,6):
            continue
        p=1 if q==6 else 3
        weights=[v.real for v in source[p-1]]
        tv=sum(abs(w) for w in weights)
        loss=tv*error
        assert loss<1
        ratios=[v/levels[-1] for v in levels[:-1]]+[arb(1)]
        if p==3:
            # 上弦界：x_j/x19 <= (x1/x19)^theta_j。
            cap=[bandwidth**(2*r.log()/ratios[0].log()) for r in ratios[:-1]]
            signed_floor=weights[-1] - sum(max(arb(0),-w.lower())*v
                                       for w,v in zip(weights[:-1],cap,strict=True))
            assert signed_floor > (arb("5.88") if q==1 else arb("2.22")) and weights[-1]>0
            print(q,p,"一阶矩除以最大位置的正下界",signed_floor.lower().str(14),flush=True)
            continue
        assert weights[-1]>0
        negative=sum(max(arb(0),-w.lower()) for w in weights[:-1])
        x18=(1-loss)/(alpha*negative)
        assert 0<x18<1
        # 向左外推末段的斜率。x19<=1 给 x_j>=x18^(log(r_j)/log(r18))。
        floor=x18**(ratios[15].log()/ratios[17].log())
        assert floor>p0
        high_error=arb("2.398e-21")+arb("7.468e-36")/floor
        max_weak=base+sum(max(arb(0),w.upper()) for w in weights[:15])*(d0-base)
        max_weak+=tv*(high_error+arb("2.398e-21")+arb("7.468e-36"))
        print(q,p,"第18层、第16层平方位置下界、全弱差上界",
              x18.lower().str(14),floor.lower().str(14),max_weak.upper().str(14),flush=True)
        assert x18>arb(".14330") and floor>arb(".01824") and max_weak<arb("-.03716")


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
