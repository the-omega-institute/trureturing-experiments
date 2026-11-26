"""Mixed-quark fit from a common 1 TeV table, uniform in both arithmetic phases."""
from pathlib import Path
import runpy

from flint import arb

root = Path(__file__).resolve().parent
central = runpy.run_path(str(root / "ch39-one-tev-quark-fit.py"))
window = runpy.run_path(str(root / "ch39-full-ckm-clock-window.py"))
metric = central["metric"]
angles = central["angle"]
radius = arb("5e-8")

def verify():
    centers = window["generators"](metric["model"])
    clocks = metric["common"]["angles"]
    varying = window["generators"](window["fit"]["model"](*(v+arb(0,radius) for v in clocks)))
    maxima = [arb(0), arb(0)]
    for i in range(metric["ncell"]):
        q = metric["qlo"]+(arb(i)+arb(1)/2)*metric["step"]
        weights = central["weights"](q)
        for j, part in enumerate((weights[:11],weights[11:])):
            error = window["perturbation"](centers,varying,q,part).upper()
            assert error.is_finite() and error > 0
            maxima[j] = max(maxima[j], error)
    dilation = (1+metric["step"]/(2*metric["qlo"]))**2
    eps = [v*dilation for v in maxima]
    targets = [metric["common"]["up"],list(reversed(metric["common"]["down"]))]
    for err,refs in zip(eps,targets):
        for ref in refs:
            rel = ((1+arb(0,".00013")) + arb(0,err)/ref)/(1+arb(0,err))
            assert abs(rel-1) < arb(".002")
    bottom = (1+arb(0,".00013"))*(1+arb(0,eps[1]))/(1+arb(0,eps[0]))
    assert abs(bottom-1) < arb(".000131")
    # Uniform absolute gaps between different singular values, in units of the largest.
    gaps = [arb(".0035"),arb(".018")]
    for (light,middle),gap in zip(targets,gaps):
        assert middle*(1-arb(".00013"))-light*(1+arb(".00013"))>gap
        assert 1-middle*(1+arb(".00013"))>gap
    dv = sum(2*err/(gap-err) for err,gap in zip(eps,gaps))
    a,b,c = [(v+arb(0,".00002")).sin() for v in angles[:3]]
    cc = (1-c*c).sqrt()
    nc = c+arb(0,dv)
    na = (a*cc+arb(0,dv))/(1-nc*nc).sqrt()
    nb = (b*cc+arb(0,dv))/(1-nc*nc).sqrt()
    angle_err = max(abs(v.asin()-ref).upper() for v,ref in zip((na,nb,nc),angles))
    min_entries = [(a*cc).lower(),(b*cc).lower(),c.lower(),
                   (((1-a*a)*(1-b*b)).sqrt()-a*b*c).lower()]
    def offset(x,y,z):
        return (x*y*z/((1-x*x)*(1-y*y)).sqrt()).asin()
    phase_err = arb(".004785")+sum((dv/v).asin() for v in min_entries)+offset(a,b,c)+offset(na,nb,nc)
    delta = angles[3]+arb(0,phase_err.upper())
    j = na*nb*nc*((1-na*na)*(1-nb*nb)).sqrt()*(1-nc*nc)*delta.sin()
    assert angle_err < arb(".00008") and phase_err < arb(".023")
    assert arb(".0000315") < j < arb(".0000331")
    print("relative matrix bounds",*(v.upper().str(10) for v in eps))
    print("CKM angle and phase errors",angle_err.str(10),phase_err.upper().str(10))
    print("J bounds",j.lower().str(12),j.upper().str(12))
    print("Haar measure",(radius/arb.pi())**2)

if __name__ == "__main__":
    verify()
