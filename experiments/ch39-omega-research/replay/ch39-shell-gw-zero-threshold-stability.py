"""Test simultaneous threshold cancellation, BFB and scalar J=0 unitarity in the 43rd window."""
from pathlib import Path
from runpy import run_path
from flint import arb

base=run_path(str(Path(__file__).with_name("ch39-shell-gw-inert-quartic-band.py")))
model,flow,q,loop=(base[k] for k in ("model","flow","q","loop"))


def check():
    _,m=flow["masses"](43)
    eps=q("12e-8")/((m[0]**2-125**2)/q("1.2e12")**2)
    b=((1+eps)/(1-eps)).sqrt()
    spread=(1/b).union(b)
    ms=[model["shifted"](v*spread,q(".36").union(q(".391"))) for v in m[:3]]
    x=q(".010")
    assert ms[1]<ms[0]*(loop*x).exp()<ms[2]
    y,_=base["matched"]["get_boundary"](ms[0].lower(),ms[0].upper())
    initial=[q(5,3)*y[2]**2,y[1]**2,y[0]**2]
    y0=sum(v*v for i,v in enumerate(y) if i not in (0,1,2,6))
    assert 10*x*y0<1
    ymax=y0/(1-10*x*y0)
    # Validated comparison bounds, allowing both active counts throughout.
    g=initial
    max_gauge=[]
    step=x/64
    for _ in range(64):
        box=[t+arb(0,".002") for t in g]
        dg=[]
        for i in range(3):
            d=model["coeff"]["db"][i]
            k=q(1).union(2)
            B=[a+k*b for a,b in zip(model["coeff"]["base"][i],model["coeff"]["extra"][i])]
            f=2*box[i]**2*(flow["beta"][i]+k*d+(sum(v*w for v,w in zip(B,box))-
                  q(0).union(max(model["coeff"]["yukawa"][i])*ymax))/loop)
            assert abs(step*f)<q(".00199")
            dg.append(f)
        max_gauge.append(box)
        g=[v+step*f for v,f in zip(g,dg)]
    g=[max(v[i].upper() for v in max_gauge) for i in range(3)]
    W,Y=g[1],q(3,5)*g[0]
    gy=sum(v*c for v,c in zip(g,flow["C"][0]))
    light=y[6].upper()
    top=y[3].lower()
    scalar=4*arb.pi()/3
    p0=q(".035")
    trace=p0
    down2=sum(y[i]**2 for i in (4,8,10))*(9*ymax*x).exp()
    gauge4=q(3,8)*(2*W**2+(W+Y)**2)
    force=q(3,2)*(3*W**2+Y**2)
    step=x/128
    hit=None
    for j in range(128):
        # First close a joint tube for the upper light quartic and the
        # largest absolute portal mass trace over this step.
        lc=light+q(".00012")
        pc=trace+q(".003")
        gap=2*(lc*scalar).sqrt()
        growth=12*(lc+scalar)+6*ymax+16*arb.pi()
        rate=growth*pc+2*pc**2+24*(pc/2+gap)**2+force
        assert step*rate<q(".003")
        crude=top*(-(gy+q("1.5")*down2)*step).exp()
        refined=top*((q("4.5")*crude**2-gy-q("1.5")*down2)*step).exp()
        lower=arb.min(top,refined)
        portal=pc**2/2+2*(pc/2+gap)**2
        rhs=24*lc**2+12*ymax*lc+gauge4-6*lower**4+2*portal
        assert step*rhs<q(".00012")
        light=(light+step*rhs).upper()
        trace=(trace+step*rate).upper()
        top=refined.lower()
        if light<0:
            hit=j+1
            break
    assert hit is not None
    print("initial |2a+b| at most",p0.str(4),"forces a negative upper quartic by GeV",
          (ms[0]*(loop*hit*step).exp()).upper().str(9),"or earlier loss of scalar J=0 bound")
    mass=p0*ms[0]**2/(16*arb.pi()**2)
    assert mass>q("1e18")
    print("one of the first two finite Higgs mass matching jumps must exceed",mass.lower().str(8))


if __name__=="__main__":
    check()
