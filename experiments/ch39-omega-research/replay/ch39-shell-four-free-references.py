"""Certify the four-shell ground index after releasing all shared length phases."""
from itertools import product
from pathlib import Path
from runpy import run_path
from flint import acb, arb

m=run_path(str(Path(__file__).with_name("ch39-shell-third-shell.py")))
penalty=m["penalty"]
layout={17:(0,1),26:(2,1),29:(2,0)}
pi=arb.pi()


def state(h):
    return {r:m["measured"](h,r) for r in m["shells"]}


def losses(obs, angles):
    return sum(m["distance"](obs[r], m["readout"](m["shells"][r][1],angles[i]*pi,angles[j]*pi))
               for r,(i,j) in layout.items())


def paired_lower(obs):
    dif=(obs[17][0]-obs[29][1],obs[17][1]-obs[26][1],obs[26][0]-obs[29][0])
    return sum(v*v for v in dif)/2


def interval_lower(obs,box):
    # Center the two equal-reference terms at their arithmetic mean to retain correlation.
    low=paired_lower(obs)
    shared=[(obs[17][0]+obs[29][1])/2,(obs[17][1]+obs[26][1])/2,(obs[26][0]+obs[29][0])/2]
    torus=[lo.union(hi)*pi for lo,hi in box]
    pairs=[x.sin_cos() for x in torus]
    for (s,_),av in zip(pairs,shared):
        low+=2*(s**2-av)**2
    for r,(i,j) in layout.items():
        s,c=pairs[i]
        t,v=pairs[j]
        mixed=1-(c*v+m["shells"][r][1]*s*t)**2
        low+=(mixed-obs[r][2])**2
    return low


def exclude_cube(obs,h,upper):
    # Common sign reversal preserves the joint readouts.
    initial=((arb(0),arb(1)/2),(arb(0),arb(1)),(arb(0),arb(1)))
    leaves=0
    def split(box):
        nonlocal leaves
        if interval_lower(obs,box)+penalty*h*h>upper:
            leaves+=1
            return
        i=max(range(3),key=lambda k:float(box[k][1]-box[k][0]))
        a,b=box[i]
        assert b-a>arb(1)/1000000
        mid=(a+b)/2
        for seg in ((a,mid),(mid,b)):
            split(box[:i]+(seg,)+box[i+1:])
    split(initial)
    return leaves


def discrete(obs):
    triple=min(losses(obs,[arb(k)/9 for k in idx]).lower() for idx in product(range(9),repeat=3))
    other=min(m["distance"](obs[867],m["reference"](867,i,j)).lower() for i,j in product(range(9),repeat=2))
    return triple+other


def main():
    h25=acb.zeta_zero(25).imag
    obs25=state(h25)
    u25=losses(obs25,[arb(k)/9 for k in (4,5,1)])+m["distance"](obs25[867],m["reference"](867,3,7))
    u9=u25+penalty*h25*h25
    assert u9<arb("0.004604")
    h719=acb.zeta_zero(719).imag
    lb719=paired_lower(state(h719))
    assert u25<arb("0.0046")<arb("0.0808")<lb719 and h25<h719
    print("all nonnegative height coefficients: 25 geometric upper",u25.str(14),"719 lower",lb719.str(14))
    # The two free R=867 phases can exactly equal their observed phases at any given height.
    h255=acb.zeta_zero(255).imag
    obs255=state(h255)
    free_upper=losses(obs255,[arb(k)/100000 for k in (83790,45394,86732)])+penalty*h255*h255
    assert free_upper<arb("0.001767")
    cap=int((u9/penalty).sqrt().upper().ceil().unique_fmpz())
    assert penalty*cap*cap>u9
    count=arb(cap).zeta_nzeros().unique_fmpz()
    assert count is not None
    maybe9=[]
    maybecont=[]
    for n,z in enumerate(acb.zeta_zeros(1,int(count)),1):
        h=z.imag
        obs=state(h)
        low=paired_lower(obs)+penalty*h*h
        if low>u9:
            continue
        maybe9.append(n)
        if n!=25:
            assert discrete(obs)+penalty*h*h>u9
        if low>free_upper:
            continue
        maybecont.append(n)
        if n!=255:
            assert n==25
            print("boxes certifying the continuous exclusion at 25",exclude_cube(obs,h,free_upper))
    assert maybe9==[25,255,309,459] and maybecont==[25,255]
    print("order-nine unique ground index 25, upper",u9.str(16))
    print("continuous unique ground index 255, upper",free_upper.str(16))
    print("uniform tail height",cap,"zeros below it",count)


if __name__=="__main__":
    main()
