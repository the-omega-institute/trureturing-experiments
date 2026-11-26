"""One-loop Abelian kinetic mixing in the window between the light neutral vector and M3."""
from pathlib import Path
from runpy import run_path

from flint import arb

root = Path(__file__).parent
flavour = run_path(str(root/"ch39-neutral-yukawa-window.py"))
study, q = flavour["study"], flavour["q"]
retuned, comp = study["retuned"], study["comp"]


def hypercharge_window(onset):
    box = None
    for index in range(3):
        for upper in (False,True):
            state = flavour["before"](index,upper)
            for _ in range(10):
                state = comp["advance"](state,(onset-q(23))/10,upper,flow=study["adjoint"]["flow"](False,False))
            value = (q("3/5")*state[0]).sqrt()
            box = value if box is None else box.union(value)
    return box


def polarization(alpha,beta,family):
    # The first row/column uses Y without its gauge coupling.
    yy,yv,vv = arb(0),arb(0),arb(0)
    for n,(d,r,sign) in zip((6,3,3,2,1),retuned["current"]["weyl_charges"].values(),strict=True):
        y = q(d)/6+q(r)/2
        for x in (1,1,-2):
            v = d*alpha+r*beta+sign*x*family
            weight = q("2/3")*n
            yy += weight*y*y
            yv += weight*y*v
            vv += weight*v*v
    # One complex Higgs doublet with Y=1/2 and 2T_R=1.
    yy,yv,vv = yy+q("1/6"),yv+beta/3,vv+2*beta**2/3
    assert (yy-q("41/6")).contains(0)
    return yy,yv,vv


def run_window(duration,gY,couplings,*,triplet=False):
    alpha,beta,family = couplings
    yy,yv,vv = polarization(*couplings)
    if triplet:
        yy += q("4/9")
        yv += q("8/3")*alpha
        vv += 16*alpha**2
    def running(t):
        p = t/(8*arb.pi()**2)
        a = 1-p*yy*gY**2
        off = -p*yv*gY
        b = 1-p*vv
        schur = b-off**2/a
        assert min(a,schur)>0
        induced = off*gY/a
        norm = 1/schur.sqrt()
        return ((alpha-induced/6)*norm,(beta-induced/2)*norm,family*norm)
    total = [arb(0)]*3
    rates = [arb(0)]*3
    charges = running(arb(0).union(duration))
    low,high,_ = flavour["charge_bounds"](*charges)
    steps = 1024
    for i in range(steps):
        t = (duration*i/steps).union(duration*(i+1)/steps)
        charges_at_t = running(t)
        increments = retuned["current"]["gauge_rates"](*charges_at_t)
        if triplet:
            extra = 8*(4*charges_at_t[0])**2/comp["loop"]**2
            increments[0] += q("4/5")*extra
            increments[2] += extra/2
        for k,(r,v) in enumerate(zip(study["adjoint"]["ratios"],increments,strict=True)):
            total[k] += v*duration/(r*steps)
            rates[k] = rates[k].union(v/r)
    return total,[v.upper() for v in rates],low,high,running(duration)


def verify():
    logs = retuned["logs"]
    h3,_,cal = map(q,logs)
    pole,old,couplings,_ = retuned["current"]["light_current"](logs)
    onset = cal+pole.log()/2
    gY = hypercharge_window(onset)
    total,rates,low,high,upper = run_window(h3-onset,gY,couplings)
    print("hypercharge at onset",gY.str(12))
    print("relative site-current increases",*((v/v0-1).str(10) for v,v0 in zip(upper,couplings)))
    print("extra gauge shifts",*((a-b).str(9) for a,b in zip(total,old)))
    params = pole,total,low,high,None
    flavour["verify_classification"](logs,params,q(".00000017"),q("4.4e-7"),rates)
    joint = run_path(str(root/"ch39-measured-neutral-yukawa.py"))
    _,_,L,z,rho,eta,c,masses = joint["fit"]["inferred_parameters"]()
    cal = study["nominal"][2]+L.log()/2
    h3,h8 = [cal+m.log() for m in masses]
    energy,observed,_,_ = joint["neutral"]["first_window"](*masses,c,z/rho,eta,1/rho)
    onset = cal+energy.log()/2
    total,rates,low,high,_ = run_window(h3-onset,hypercharge_window(onset),observed)
    assert max(rates)<q(".001")
    joint["exclude"](onset,h3,h8,low,high,total[1],total[0]-total[1])


if __name__ == "__main__":
    verify()
