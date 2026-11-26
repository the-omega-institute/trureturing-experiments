"""在三段时钟中让中心占比自由变化，证明严格误差预算下的最小总拍数。"""
from fractions import Fraction
from itertools import permutations
from pathlib import Path
from runpy import run_path

from flint import acb, arb

root=Path(__file__).parent
global_clock=run_path(str(root/"ch39-global-three-zero-clock.py"))
minimum=run_path(str(root/"ch39-ordered-three-zero-minimum.py"))
fixed,edges,ordered,support,study,box=(global_clock[k] for k in ("fixed","edges","ordered","support","study","box"))
bound=725185


def fractions_in_interval(lo,hi):
    # 先找开区间内的最简分数，再沿固定阶的法里数列分别向左右枚举。
    _,_,p,q,_,_=fixed["least_clock"](lo,hi)
    if q>bound:
        return
    middle=Fraction(p,q)
    yield middle
    for direction,edge in ((-1,lo),(1,hi)):
        f=middle
        while True:
            p,q=f.numerator,f.denominator
            d=(-direction*pow(p,-1,q))%q
            d+=((bound-d)//q)*q
            c=(p*d+direction)//q
            assert c*q-p*d==direction and q+d>bound and 0<d<=bound
            f=Fraction(c,d)
            diff=direction*(arb(c)/d-edge)
            assert 0 not in diff
            if diff>0:
                break
            yield f


def verify():
    cost,_,scale,*_=support["prev"]["setup"]()
    r0=arb(1)/3
    cache={}
    def atom(n):
        if n not in cache:
            height=acb.zeta_zero(n).imag
            x=study["center"]/(height**2+study["offset"])-1 if n!=11431 else arb(0)
            c=cost(x)[0] if n!=11431 else arb(0)
            cache[n]=x,c,study["native"]["plane_response"](height)
        return cache[n]
    corner=[v+s*arb(2)**-28 for v,s in zip(box["centers"],(1,-1,1,-1))]
    data=list(map(atom,(11431,11432,11415,11430,11372,11304)))
    ws,_=edges["candidate"]([x+scale*c/2 for x,c,t in data],
                            [box["flavor_data"](t,corner)[0] for x,c,t in data],r0)
    cost_limit=support["dot"](ws,[c for x,c,t in data])
    e=scale*cost_limit/2
    def lower_cost(l,h):
        x,c,_=atom(l)
        y,d,_=atom(h)
        return r0*d+((-e).union(e)-r0*y)/(x-y)*(c-d)
    assert r0*atom(11430)[0]>e and r0*atom(11432)[0]<-e
    assert lower_cost(11380,11432)>cost_limit and lower_cost(11430,11510)>cost_limit
    cube=[v+arb(0,arb(2)**-28) for v in box["centers"]]
    tolerances=[arb(".000055")]*3+[arb(".006")]
    winners=[]
    cases=0
    for l in range(11381,11431):
        a,c,tl=atom(l)
        for h in range(11432,11510):
            if lower_cost(l,h)>cost_limit:
                break
            b,d,th=atom(h)
            lowers,uppers=[],[]
            for x,y in ((a,b),(a+scale*c,b+scale*d)):
                lowers.append((-e/r0-y)/(x-y))
                uppers.append((e/r0-y)/(x-y))
            lo,hi=max(lowers),min(uppers)
            assert 0 not in lo-hi
            if lo>hi:
                continue
            assert 0<lo<hi<1
            ns=(11431,l,h)
            mats=(atom(11431)[2],tl,th)
            modes=[minimum["spectral_modes"](mat,cube) for mat in mats]
            for frac in fractions_in_interval(lo,hi):
                k,m=frac.numerator,frac.denominator
                f=arb(k)/m
                g=b+f*(a-b)
                z=g+scale*(d+f*(c-d))
                per_external=max(abs(g),abs(z))/e
                for t in range(1,bound//m+1):
                    ext=t*m
                    cutoff=(ext*per_external).floor().unique_fmpz()
                    assert cutoff is not None
                    shortest=max(ext+1,int(cutoff)+1)
                    longest=min(3*ext,bound)
                    if shortest>longest:
                        continue
                    cases+=1
                    start,end=arb(ext)/longest,arb(ext)/shortest
                    for order in permutations(range(3)):
                        if (l,h,t*k,ext,shortest,order)==(11410,11432,10985,338487,bound,(2,0,1)):
                            assert shortest==longest==bound
                            weights=(1-start,start*f,start*(1-f))
                            steps={n:ordered["evolution"](mat,w,cube) for n,mat,w in zip(ns,mats,weights)}
                            obs=ordered["readings"](steps,tuple(ns[i] for i in order))
                            assert all(abs(o-v)<d for o,v,d in zip(obs,ordered["flavor"]["shell"]["centers"],tolerances))
                            winners.append(tuple(ns[i] for i in order))
                            print("Winning shifts X,Z",(start*g).str(19),(start*z).str(19),"CKM",*[o.str(10) for o in obs],flush=True)
                        else:
                            minimum["exclusion"](order,ns,modes,(f,1-f),end,tolerances,start)
    assert winners==[(11432,11431,11410)]
    print("Surviving integer count cases before flavor",cases,"unique total",bound)


if __name__=="__main__":
    verify()
