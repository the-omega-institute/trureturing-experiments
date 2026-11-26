"""All 53 finite threshold jumps and one-loop mass transport in the inert reference model."""
from pathlib import Path
from runpy import run_path
from flint import arb, arb_mat

base=run_path(str(Path(__file__).with_name("ch39-shell-gw-full-inert-potential.py")))
q,loop=base["q"],base["loop"]
unit=q("1e24")


def portal(quart,i,j):
    ij=base["pair"](i,j)
    return 2*quart["a",ij]+quart["b",ij]


def threshold(k,classes,quart,m,values,derivative):
    group=1 if k<=2 else 2
    if len(values)<len(classes):
        values.append(q(0))
        derivative.append(q(0))
    mass=m**2/unit
    for i,n in classes.items():
        n-=int(i==group)
        if n:
            values[i]+=n*portal(quart,i,group)*mass/loop
    values[group]+=mass
    assert min(v.lower() for v in values)>0


def transport(classes,quart,gs,ym,dt,values,derivative):
    n=len(classes)
    G=arb_mat(n,n)
    damp=q(3,2)*(3*gs[1]+q(3,5)*gs[0])
    for i,num in classes.items():
        G[i,i]=12*quart["l",i]-damp+(q(0).union(6*ym) if i==0 else 0)
        if num>1:
            G[i,i]+=2*(num-1)*portal(quart,i,i)
        for j in classes:
            if i!=j:
                G[i,j]=2*num*portal(quart,i,j)
                assert G[i,j]>0
    lo=arb_mat(n,n,[v.lower() for v in G.entries()])
    hi=arb_mat(n,n,[v.upper() for v in G.entries()])
    plo,phi=(dt*lo).exp(),(dt*hi).exp()
    # Entrywise order is valid for Metzler systems with nonnegative inputs.
    for vector in (values,derivative):
        lower=plo*arb_mat([[arb.max(v.lower(),q(0))] for v in vector])
        upper=phi*arb_mat([[v.upper()] for v in vector])
        vector[:]=[arb.max(a.lower(),q(0)).union(b.upper()) for a,b in zip(lower.entries(),upper.entries())]


def check():
    events=[]
    base["check"](mass_path=events)
    low=q("-1e6").union(q("1e6"))
    values,derivative=[low/unit],[q(1)]
    light_after_two=None
    for op,*args in events:
        if op=="threshold":
            threshold(*args,values,derivative)
            if args[0]==2:
                light_after_two=values[0]*unit
        else:
            transport(*args,values,derivative)
    light=values[0]*unit
    slope=derivative[0]
    precision=slope*(q(125)**2/2)/light
    print("light quadratic coefficient GeV² above second and at Q",
          [(float(v.lower()),float(v.upper())) for v in (light_after_two,light)],
          "response of high light to low light",float(slope.lower()),float(slope.upper()),
          "weak-size relative change at Q below",precision.upper().str(8))
    assert light>q("1e23") and slope>0 and precision<q("1e-19")


if __name__=="__main__":
    check()
