"""总带电中性质量差为零时，十六层固定重谱的中性残差上界。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import arb, ctx


def verify(root):
    ctx.prec = 192
    energy, portal = map(arb, ("1.176397e14", ".42995401872"))
    system = run_path(str(root/"ch39-six-broken-conditions.py"))["broken_system"](root, energy, portal)
    weights = system["weights_at"](portal)[0]
    data = system["original"]
    sites, p0, p1 = data["sites"], data["p0"], data["p1"]
    centers = [z+s*portal*system["rho"] for z,s in zip(system["masses"],system["signs"],strict=True)]
    masses = [z.sqrt() for z in centers]
    assert len(masses)==16 and all(z>0 for z in masses)
    assert all(weights[j]>0 for j in (0,2,4)) and all(weights[j]<0 for j in (1,3))
    assert sites[4] > sites[2]
    run_path(str(root/"ch39-neutral-residual-all-cutoffs.py"))["verify"]()
    light = arb("-.074256577194663365546").union(arb("-.074256577194663365544"))
    light += arb(0,arb("4.907e-20"))
    gauge = run_path(str(root/"ch39-all-field-signed-crossing.py"))
    rule = gauge["angular_rule"](24,57)
    log_ratio = (p1/p0).log()
    multiplier = arb("12.6")
    for value, target in ((".3941","-.00038308"),(".3955","-.00003377"),(".3956","-.000008787")):
        limit = arb(value)
        cap = max(z.upper() for z in centers)*(1+limit)**2+p1/4
        ratio = cap/sites[2]
        floor = (weights[2]+weights[4])*(1+ratio)*(-ratio).exp()+weights[1]+weights[3]
        assert floor > arb(".8")
        full_residual = light
        penalty = arb(0)
        for center, mass in zip(centers,masses,strict=True):
            def profile(t):
                return t/(center*t*t+p1/4)**2
            endpoints = [profile(t).lower() for t in (1-limit,arb(1),1+limit)]
            slope = (arb(".8")*center*(p1-p0)/(240*arb.pi()) *
                     (min(endpoints[:2])+min(endpoints[1:])))
            penalty += 2*limit*max(arb(0),(multiplier*mass-slope).upper())
            for cutoff, weight in zip(sites,weights,strict=True):
                left, right = [gauge["angular_exposure"](
                    1,center*(1+d)**2,cutoff,p0,p1,log_ratio,rule) for d in (-limit,limit)]
                full_residual += weight*(left-right)/(48*arb.pi())
        bound = full_residual-multiplier*limit*sum(masses)+penalty
        assert bound < arb(target)
        print("幅度上限及总质量差为零时的残差上界", value, bound.upper().str(12))
        if value==".3941":
            required = -bound/multiplier
            assert required > arb(".000030403")
            assert 2*arb("1.22089e19")*required > arb("7.42e14")


if __name__=="__main__":
    verify(Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent)
