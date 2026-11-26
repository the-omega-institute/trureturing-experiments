"""One-loop threshold test for electroweak doublets on the eight lowest shell modes."""
from pathlib import Path
from runpy import run_path
import numpy as np
from flint import arb, arb_mat

model=run_path(str(Path(__file__).with_name("ch39-shell-gw-195-flavor.py")))
rows=model["bands"]["m"]["rows"]
count=model["count"]
mz=arb("91.1876")


def crossing():
    em,s=arb("127.930"),arb(".23129")
    a1,a2,a3=3*(1-s)*em/5,s*em,1/arb(".1180")
    b1,b2,b3=arb(41)/10,-arb(19)/6,arb(-7)
    h1,h2=arb(1)/10,arb(1)/6
    coeff=arb_mat([[b1-b2,h1-h2],[b2-b3,h2]])
    L,E=coeff.solve(2*arb.pi()*arb_mat([[a1-a2],[a2-a3]])).entries()
    assert arb("6.4451")<E/L<arb("6.4453")
    Q=mz*L.exp()
    assert arb("5.0144e13")<Q<arb("5.0145e13")
    return Q,E,a3+7*L/(2*arb.pi())


def lowest_eight(index):
    if index==43:
        sing,family,_,_=model["m"]["transported_spectra"]()
    else:
        sing,family,_=model["scalar_models"]()
    s,f=rows(sing),rows(family)
    proposals=sorted(np.r_[np.linalg.eigvalsh(model["array"](sing)),
                            np.linalg.eigvalsh(model["array"](family))])[:8]
    levels=[]
    for i,guess in enumerate(proposals):
        shift=arb("1e-11")
        lo,hi=arb(float(guess))-shift,arb(float(guess))+shift
        for c,k in ((lo,i),(hi,i+1)):
            assert count(s,c)+count(f,c)==k
        levels.append(lo.union(hi))
    return [e-levels[0] for e in levels[1:]]


def check(index,Q,E):
    levels=lowest_eight(index)
    def masses(scale,shift):
        return [(shift**2+scale**2*d).sqrt() for d in levels]
    # One common additive squared-mass input, given in GeV^2.
    shift=arb(125)
    lo,hi=map(arb,("1471.22","1471.25") if index==43 else ("1297.88","1297.92"))
    assert sum((Q/m).log() for m in masses(lo,shift))>E
    assert sum((Q/m).log() for m in masses(hi,shift))<E
    actual=masses(lo.union(hi),shift)
    assert actual[0]>mz and actual[-1]<Q
    zero_scale=Q*(-E/7).exp()*(-sum(z.log() for z in levels)/14).exp()
    light=zero_scale*levels[0].sqrt()
    assert light< arb(84) if index==43 else light>arb(162)
    print("zero",index,"mass scale",lo.union(hi).str(9),"GeV, lightest and heaviest added masses",
          actual[0].str(8),actual[-1].str(8),"unshifted lightest",light.str(9))


if __name__=="__main__":
    Q,E,alpha=crossing()
    print("one-loop Q in GeV",Q.str(12),"total extra log exposure",E.str(12),"inverse coupling",alpha.str(10))
    for n in (43,195):
        check(n,Q,E)
