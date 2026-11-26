"""A signed lower bound for the pair of near-degenerate heavy mass thresholds."""
from pathlib import Path
from runpy import run_path
from flint import arb

base=run_path(str(Path(__file__).with_name("ch39-shell-gw-inert-quartic-band.py")))
model,flow,q,loop=(base[k] for k in ("model","flow","q","loop"))


def check():
    _,m=flow["masses"](43)
    eps=q("12e-8")/((m[0]**2-125**2)/q("1.2e12")**2)
    b=((1+eps)/(1-eps)).sqrt()
    masses=[v*(1/b).union(b) for v in m[:3]]
    s0,s1=q(".36"),q(".391")
    def logratio(i,j):
        def at(s):
            return (model["shifted"](masses[i],s)/model["shifted"](masses[j],s)).log()
        return at(s0).union(at(s1))
    split=logratio(1,0)/loop
    third=logratio(2,0)/loop
    ratio=(2*loop*split).exp()
    first=model["shifted"](masses[0],s0.union(s1))
    y,_=base["matched"]["get_boundary"](first.lower(),first.upper())
    top0=y[3].lower()
    total0=sum(v*v for i,v in enumerate(y) if i not in (0,1,2,6))
    end=q(".012")
    assert 0<split<end<third
    total=total0/(1-10*end*total0)
    gu=[q(5,3)*y[2]**2,y[1]**2,y[0]**2]
    # Gauge bounds to the third threshold; permit either one or two extra scalars.
    h=third.upper()/128
    allg=[]
    for _ in range(128):
        box=[t+arb(0,".003") for t in gu]
        dg=[]
        for i in range(3):
            k=q(1).union(2)
            row=[v+k*a for v,a in zip(model["coeff"]["base"][i],model["coeff"]["extra"][i])]
            max_t=total0/(1-10*third.upper()*total0)
            f=2*box[i]**2*(flow["beta"][i]+k*model["coeff"]["db"][i]+
                (sum(a*g for a,g in zip(row,box))-q(0).union(max(model["coeff"]["yukawa"][i])*max_t))/loop)
            assert abs(h*f)<q(".00299")
            dg.append(f)
        allg.append(box)
        gu=[a+h*f for a,f in zip(gu,dg)]
    g=[max(row[i].upper() for row in allg) for i in range(3)]
    W,Y=g[1],q(3,5)*g[0]
    assert 3*(3*W+Y)<4
    G=q(3,8)*(2*W**2+(W+Y)**2)
    F=q(3,2)*(3*W**2+Y**2)
    gy=sum(v*c for v,c in zip(g,flow["C"][0]))
    down2=sum(y[i]**2 for i in (4,8,10))*(9*total*end).exp()
    def cap(x):
        # l_j' >= 24*l_j²-4*l_j and l_j <= 4*pi/3 at the third threshold.
        return 4/(24-(24-4/(4*arb.pi()/3))*(-4*(third.lower()-x)).exp())
    t=split.upper()
    maxL=cap(t).upper()
    l0=y[6].upper()
    assert l0<q(".00110")
    # Assume p1_at_entry+(M2/M1)^2*p2_at_entry <= allowance.
    allowance=q(".020")
    # First bound l0 at M2 without using the assumed sum. Scalar unitarity
    # bounds each portal coefficient: p1 is a trace-state matrix element;
    # a and a+b are also controlled by boundedness and p1.
    def before(pmax,light,check_tube):
        gap=2*(light*maxL).sqrt()
        d=pmax/2+gap
        pl=pmax**2/2+2*d*d
        beta=24*light**2+12*total*light+G+pl
        if check_tube:
            assert l0+t*beta<light
        return l0+t*beta
    preliminary=before(8*arb.pi(),q(".01"),True)
    p1max=allowance+ratio*4*(preliminary*maxL).sqrt()
    p1min=-4*(l0*maxL).sqrt()
    pbound=max(abs(p1min).upper(),p1max.upper())+q(".001")
    # Bootstrap the tiny first layer using its now sharper initial portal bound.
    light=q(".00112")
    gap=2*(light*maxL).sqrt()
    prate=(12*(light+maxL)+6*total+4)*pbound+2*pbound**2+24*(pbound/2+gap)**2+F
    assert t*prate<q(".001")
    light=before(pbound,light,True).upper()
    # The sum at M2 includes the drift of p1 and the small mass-ratio term.
    c=allowance+t*prate+(ratio-1)*4*(light*maxL).sqrt()
    top=top0*(-(gy+q("1.5")*down2)*t).exp()
    step=(end-t)/192
    hit=None
    for j in range(192):
        future=cap(t+(j+1)*step).upper()
        lc=light+q(".00005")
        qc=c+q(".002")
        gap=2*(lc*future).sqrt()
        # On the positive upper barrier for the trace sum, the negative
        # trace cannot be below -2*gap. These maxima use both signs.
        sq=qc**2+4*gap*qc+8*gap**2
        rate=(12*lc+6*total+16*arb.pi())*qc+12*future*(qc+2*gap)+8*sq+24*gap*qc+48*gap**2+2*F
        assert step*rate<q(".002")
        portal=qc**2+6*gap*qc+12*gap**2
        crude=top*(-(gy+q("1.5")*down2)*step).exp()
        nexttop=top*((q("4.5")*crude**2-gy-q("1.5")*down2)*step).exp()
        rhs=24*lc**2+12*total*lc+G-6*arb.min(top,nexttop)**4+portal
        assert step*rhs<q(".00005")
        light=(light+step*rhs).upper()
        c=(c+step*rate).upper()
        top=nexttop.lower()
        if light<0:
            hit=j+1
            break
    assert hit is not None
    by=(first*(loop*(t+hit*step)).exp()).upper()
    size=allowance*first**2/loop
    print("signed pair trace in first mass units must exceed",allowance.str(4),
          "negative Higgs finite steps have signed total magnitude above",size.lower().str(8),
          "contradiction by GeV",by.str(9))


if __name__=="__main__":
    check()
