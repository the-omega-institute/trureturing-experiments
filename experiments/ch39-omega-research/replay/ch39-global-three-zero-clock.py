"""固定中心占三分之二、切线误差平分时，在全部三编号单段时间表上找最短分辨时钟。"""
from itertools import permutations
from pathlib import Path
from runpy import run_path

from flint import acb, arb

fixed=run_path(str(Path(__file__).with_name("ch39-three-zero-integer-clock.py")))
support,ordered,study,box,edges=(fixed[key] for key in ("support","ordered","study","box","edges"))


def verify():
    cost,strength,scale,*_=support["prev"]["setup"]()
    center=11431
    r=arb(1)/3
    cache={}
    def atom(n):
        if n not in cache:
            height=acb.zeta_zero(n).imag
            x=arb(0) if n==center else study["center"]/(height**2+study["offset"])-1
            c=arb(0) if n==center else cost(x)[0]
            t=study["native"]["plane_response"](height)
            cache[n]=(x,c,t)
        return cache[n]
    corner=[v+s*arb(2)**-28 for v,s in zip(box["centers"],(1,-1,1,-1))]
    data=list(map(atom,(11431,11432,11415,11430,11372,11304)))
    ws,_=edges["candidate"]([x+scale*c/2 for x,c,t in data],
                            [box["flavor_data"](t,corner)[0] for x,c,t in data],r)
    cost_limit=support["dot"](ws,[c for x,c,t in data])
    limit=scale*cost_limit/2
    assert r*atom(center-1)[0]>limit and r*atom(center+1)[0]<-limit
    def min_chord(l,rgt):
        x,c,_=atom(l)
        y,d,_=atom(rgt)
        gauge=(-limit).union(limit)
        wl=(gauge-r*y)/(x-y)
        return r*d+wl*(c-d)
    low,high=11380,11510
    assert min_chord(low,center+1)>cost_limit and min_chord(center-1,high)>cost_limit
    cube=[c+arb(0,arb(2)**-28) for c in box["centers"]]
    tolerances=[arb(".000055")]*3+[arb(".006")]
    feasible=[]
    for l in range(low+1,center):
        x,c,tl=atom(l)
        for h in range(center+1,high):
            if min_chord(l,h)>cost_limit:
                break
            y,d,th=atom(h)
            lower=[]
            upper=[]
            for left,right in ((x,y),(x+scale*c,y+scale*d)):
                assert left>right
                lower.append((-limit/r-right)/(left-right))
                upper.append((limit/r-right)/(left-right))
            lo,hi=max(lower),min(upper)
            assert 0 not in lo-hi
            if lo>hi:
                continue
            assert 0<lo<hi<1
            f=lo.union(hi)
            steps={l:ordered["evolution"](tl,r*f,cube),
                   center:ordered["evolution"](atom(center)[2],1-r,cube),
                   h:ordered["evolution"](th,r*(1-f),cube)}
            for order in permutations((l,center,h)):
                obs=ordered["readings"](steps,order)
                errors=[abs(v-target) for v,target in zip(obs,ordered["flavor"]["shell"]["centers"])]
                if all(e<cut for e,cut in zip(errors,tolerances)):
                    a,b,p,q,c0,d0=fixed["least_clock"](lo,hi)
                    bits,k,den=fixed["binary_clock"](lo,hi)
                    feasible.append((q,p,order,bits,k))
                else:
                    assert any(e>cut for e,cut in zip(errors,tolerances)),(order,errors)
    for q,p,order,bits,k in sorted(feasible):
        print("Left-to-right",order,"shortest total",3*q,"external left",p,"denominator",q,
              "binary bits",bits,"left ticks",k)
    expected=[(338487,10985,11410,35),(478549,11248,11402,38),(901378,20813,11401,39),
              (1676772,36683,11399,41),(1884001,42017,11400,39)]
    assert [(q,p,order[2],bits) for q,p,order,bits,k in sorted(feasible)]==expected
    assert all(order[:2]==(11432,11431) for q,p,order,bits,k in feasible)
    print("Number of admissible triples and orders",len(feasible))


if __name__=="__main__":
    verify()
