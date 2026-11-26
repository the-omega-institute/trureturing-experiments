"""Full shell-doublet thresholds and the cost of inducing the top hierarchy in that basis."""
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import arb

eight=run_path(str(Path(__file__).with_name("ch39-shell-gw-eight-gauge-threshold.py")))
common=run_path(str(Path(__file__).with_name("ch39-shell-gw-195-common-field.py")))
model=eight["model"]
rows,count=eight["rows"],eight["count"]


def all_levels(sing,family):
    blocks=[rows(z) for z in (sing,family)]
    trials=sorted(np.concatenate([np.linalg.eigvalsh(model["array"](z)) for z in (sing,family)]))
    # Keep entire numerical clusters and verify their joint multiplicity by inertia.
    groups=[]
    for val in trials:
        if not groups or val-groups[-1][-1]>1e-8:
            groups.append([])
        groups[-1].append(val)
    values=[]
    for group in groups:
        lo=arb(float(group[0]))-arb("1e-10")
        hi=arb(float(group[-1]))+arb("1e-10")
        assert sum(count(b,lo) for b in blocks)==len(values)
        assert sum(count(b,hi) for b in blocks)==len(values)+len(group)
        values.extend([lo.union(hi)]*len(group))
    assert len(values)==54 and values[1]>values[0]
    return [v-values[0] for v in values[1:]]


def common_bound(n):
    if n==43:
        kin=run_path(str(Path(__file__).with_name("ch39-shell-gw-native-kinetic-cost.py")))
        K,F=kin["scalar_kinetic"](),kin["m"]["F"]
        x,lower=arb(340),arb("2.049")
    else:
        K,F=common["operators"]()
        x,lower=arb(350),arb("2.0528")
    t=common["native"]["relative_phase"](arb(676))
    Y=common["conjugate"](F,t)
    G,D=model["adj"](F)*F,model["adj"](Y)*Y
    up,down=common["old"]["u2"],common["old"]["d2"]
    assert count(rows(K+23000*(D-down/up*G)-x*G),lower-x*up/(4*arb.pi())**2)==0
    return lower


def check(n,Q,E):
    if n==43:
        single,family,_,_=model["m"]["transported_spectra"]()
        B=model["bands"]["orb"]["orthonormalize"](model["bands"]["images"](2))[0]
        F=common["old"]["F"]
        bracket=("1.11863e12","1.11867e12")
        gap=arb("5.15e10")
    else:
        single,family,F=model["scalar_models"]()
        B=model["image_basis"](F,model["sources"](single,family,2))[0]
        bracket=("1.08054e12","1.08058e12")
        gap=arb("1.25e11")
    levels=all_levels(single,family)
    def masses(scale):
        return [(125**2+scale**2*z).sqrt() for z in levels]
    lo,hi=map(arb,bracket)
    assert sum((Q/m).log() for m in masses(lo))>E
    assert sum((Q/m).log() for m in masses(hi))<E
    assert masses(lo)[0]>gap and masses(hi)[-1]<Q
    adj=model["adj"]
    assert count(rows(B*adj(B)),arb(".578")**2)==3
    assert count(rows(F*adj(F)),arb(".52")**2)==9
    top,charm,light=(arb(s) for s in (".8616",".00311",".00000615"))
    miss=(top**2+charm**2+light**2-(arb(".578")*(top+charm+light))**2).sqrt()
    tail=miss/(4*arb.pi()*arb(".52"))
    assert tail>arb(".107")
    ratio=tail*(gap/arb("246.1"))**2
    assert ratio > arb("4.6e15" if n==43 else "2.7e16")
    assert count(rows(single),arb(".014"))==1 and count(rows(family),arb(".014"))==0
    gap_kinetic=common_bound(n)-arb(".014")
    gradient=lo**2*gap_kinetic/arb("246.1")**2
    assert gradient>arb("4.20e19" if n==43 else "3.92e19")
    print("joint gradient / VEV^3 above",gradient.lower().str(7),
          "required isospectral change / M² above",gap_kinetic.str(6))
    print("zero",n,"full-spectrum mass scale in",bracket,"first heavy mass above",gap.str(5),
          "VEV fraction needed above",tail.lower().str(7),"gradient / VEV^3 above",ratio.lower().str(7))


if __name__=="__main__":
    Q,E,_=eight["crossing"]()
    for index in (43,195):
        check(index,Q,E)
