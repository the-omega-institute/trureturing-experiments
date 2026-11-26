"""Enclose open arithmetic windows on which all CP phases coexist at fixed coupling."""
from pathlib import Path
from runpy import run_path
from flint import acb, acb_mat, arb, ctx

cp=run_path(str(Path(__file__).with_name("ch39-shell-gw-fixed-strength-cp.py")))
model=cp["model"]
sel=model["selected"]
rows=cp["rows"]; count=model["count"]


def window(index,width):
    ctx.prec=512
    h=acb.zeta_zero(index).imag
    p,q=[(a-model["angles"][0])*h for a in model["angles"][1:]]
    P,Q=p+arb(0,width),q+arb(0,width)
    width_error=[]
    for sign in (1,-1,0):
        _,center=sel["action_kernel"](sign,p,q)
        _,box=sel["action_kernel"](sign,P,Q)
        gap=arb(".02" if index==43 else ".007")
        assert count(center,-gap)==count(center,gap)==18-3*sign
        K=acb_mat(center)
        scale=arb(8)
        assert max(sum(abs(z) for z in row) for row in center)<scale
        # Differences have entries containing zero; enclose the Frobenius
        # norm by magnitude summation if a complex-box dot product straddles 0.
        delta=sum(abs(z).abs_upper()**2 for z in (acb_mat(box)-K).entries()).sqrt()
        assert delta<gap
        move=delta/(2*(gap-delta))
        x=K/scale
        b=gap/scale
        unit=model["eye"](36)
        while 1-b>arb("1e-16"):
            x=x*(3*unit-x*x)/2
            b=b*(3-b*b)/2
        err=arb(0,(1-b).upper()/2)
        sep=arb(".15" if sign else ".06")
        for i,s in ((0,1),(18,-1)):
            block=acb_mat([[(unit[j,k]-s*x[j,k])/2+acb(err,err)
                            for k in range(i,i+18)] for j in range(i,i+18)])
            assert count(rows(block),sep)==(3 if s==sign else 0)
        assert move<sep
        if sign:
            width_error.append(arb(2).sqrt()*move/sep)
    return arb(2).sqrt()*sum(width_error)


def gauge_window(index,width,Q,E):
    gauge=run_path(str(Path(__file__).with_name("ch39-shell-gw-full-gauge-flavor.py")))
    if index==43:
        sing,family,_,_=model["m"]["transported_spectra"]()
    else:
        sing,family,_=model["scalar_models"]()
    gaps=gauge["all_levels"](sing,family)
    assert gaps[0]>4*width
    enclosed=[gap+arb(0,4*width) for gap in gaps]
    def masses(scale):
        return [(125**2+scale**2*gap).sqrt() for gap in enclosed]
    low=masses(arb("1e12"))
    high=masses(arb("1.2e12"))
    assert min(low)>arb("91.1876") and max(high)<Q
    assert sum((Q/a).log() for a in low)>E
    assert sum((Q/a).log() for a in high)<E


if __name__=="__main__":
    widths=[arb("3e-8"),arb("1e-8")]
    gauge=run_path(str(Path(__file__).with_name("ch39-shell-gw-eight-gauge-threshold.py")))
    Q,E,_=gauge["crossing"]()
    for n,width in zip((43,195),widths):
        difference=window(n,width)
        cp["check"](n,difference)
        gauge_window(n,width,Q,E)
        print("zero",n,"joint output difference below",difference.upper().str(8))
    h43,h195=(acb.zeta_zero(n).imag for n in (43,195))
    dp=(model["angles"][1]-model["angles"][0])*(h195-h43)
    assert abs((dp/2).sin())>sum(widths)/2
    density=sum((w/arb.pi())**2 for w in widths)
    print("lower relative density above",density.lower().str(9))
