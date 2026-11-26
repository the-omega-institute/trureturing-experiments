"""由四项 CKM 的任意一项区分五种首段编号，并恢复连续校准与时长。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat

root=Path(__file__).parent
identification=run_path(str(root/"ch39-three-clock-flavor-identification.py"))
fixed=run_path(str(root/"ch39-three-zero-integer-clock.py"))
support,ordered,study,box,edges=(fixed[k] for k in ("support","ordered","study","box","edges"))
labels=(11410,11402,11401,11399,11400)


def verify():
    cost,_,scale,*_=support["prev"]["setup"]()
    centers=box["centers"]
    delta=arb(2)**-28
    cube=[c+arb(0,delta) for c in centers]
    r=arb(1)/3
    cache={}
    def atom(n):
        if n not in cache:
            height=acb.zeta_zero(n).imag
            x=arb(0) if n==11431 else study["center"]/(height**2+study["offset"])-1
            c=arb(0) if n==11431 else cost(x)[0]
            mat=study["native"]["plane_response"](height)
            cache[n]=(x,c,mat)
        return cache[n]
    corner=[v+s*delta for v,s in zip(centers,(1,-1,1,-1))]
    data=list(map(atom,(11431,11432,11415,11430,11372,11304)))
    w,_=edges["candidate"]([x+scale*c/2 for x,c,t in data],
                            [box["flavor_data"](t,corner)[0] for x,c,t in data],r)
    e=scale*support["dot"](w,[c for x,c,t in data])/2
    flavor=ordered["flavor"]
    ident=arb_mat([[int(i==j) for j in range(4)] for i in range(4)])
    envelopes={}
    max_eta=arb(0)
    slopes=(r,arb(0),-r)
    for n in labels:
        x,c,tl=atom(n)
        y,d,tr=atom(11432)
        lows,highs=[],[]
        for left,right in ((x,y),(x+scale*c,y+scale*d)):
            lows.append((-e/r-right)/(left-right))
            highs.append((e/r-right)/(left-right))
        lo,hi=max(lows),min(highs)
        assert lo<hi
        share=lo.union(hi)
        terms=[(v,flavor["arcs"][32],flavor["arcs"][152],flavor["identity"]) for v in (tl,atom(11431)[2],tr)]
        def inputs(f):
            return list(zip(terms,(r*f,1-r,r*(1-f))))
        vals,jac=identification["observations"](inputs(share),cube,slopes)
        a,b,u,j=vals
        one=1-u
        central=identification["observations"](inputs(share.mid()),centers)[1]
        pre=arb_mat([[central[k,i].mid() for i in range(4)] for k in range(4)]).inv()
        rem=ident-pre*arb_mat([[jac[k,i] for i in range(4)] for k in range(4)])
        eta=max(sum(abs(rem[k,i]) for i in range(4)).upper() for k in range(4))
        max_eta=max(max_eta,eta)
        assert eta<arb(".002")
        signs=[[] for _ in range(4)]
        for k in range(5):
            da,db,du,dj=(jac[i,k] for i in range(4))
            phase=2*dj/j-da/a-db/b-du/u+(da+du)/(one-a)+(db+du)/(one-b)-2*du/one
            changes=[one*da+a*du,one*db+b*du,du,phase]
            for seq,v in zip(signs,changes):
                assert 0 not in v,(n,k,v.str(8))
                seq.append(1 if v>0 else -1)
        def observe(f,params):
            ns=(n,11431,11432)
            steps={k:ordered["evolution"](row[0],duration,params) for k,(row,duration) in zip(ns,inputs(f))}
            return ordered["readings"](steps,tuple(reversed(ns)))
        ends=[]
        for k,sgns in enumerate(signs):
            limits=[]
            for side in (-1,1):
                params=[c+side*t*delta for c,t in zip(centers,sgns)]
                f=hi if side*sgns[4]>0 else lo
                limits.append(observe(f,params)[k])
            assert limits[0]<limits[1]
            ends.append(limits)
        envelopes[n]=ends
    print("Uniform inverse bound on all five labels",max_eta.str(12))
    for k in range(4):
        gaps=[]
        for m,n in combinations(labels,2):
            left,right=envelopes[m][k],envelopes[n][k]
            gap=max(left[0]-right[1],right[0]-left[1])/2
            assert gap>0
            gaps.append((gap,m,n))
        best=min(gaps,key=lambda row:float(row[0].mid()))
        assert all(best[0]<val for val,a,b in gaps if (a,b)!=best[1:])
        print("Column",k,"critical single-reading error",best[0].str(17),"pair",best[1:])


if __name__=="__main__":
    verify()
