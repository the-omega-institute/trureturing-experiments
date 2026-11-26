"""A cyclic orbit of three existing clock edges creates robust degree three at zero 211."""
from pathlib import Path
import runpy
from flint import acb,arb,acb_mat

work=Path(__file__).resolve().parent
fit=runpy.run_path(str(work/"ch39-hierarchy-in-zero-flux.py"))
creation=runpy.run_path(str(work/"ch39-source-created-flux.py"))
src=fit["src"]
s,W,J,block=src["shell"],src["W"],src["J"],src["block"]
height=acb.zeta_zero(211).imag
p,q=[s[k]*height for k in ("a","b")]
orbit=[(0,7),(3,8),(4,6)]


def source_data():
    initial,_=fit["model"](p,q)
    G=block(initial,J,J).inv()
    proj=block(initial,W,J)*G
    k=block(initial,W,W)-proj*block(initial,J,W)
    cols=[]
    for i,j in orbit:
        v=-block(proj,range(7),(J.index(i),))
        v[W.index(j),0]+=initial[i,j].conjugate()
        cols.append(v)
    return k,G,cols


def check_pair(a,b):
    k,g,v=source_data()
    coeff=[]
    for mask in ((),(a,),(b,),(a,b)):
        if not mask:
            coeff.append(k)
            continue
        ids=[J.index(orbit[i][0]) for i in mask]
        gg=block(g,ids,ids)
        vv=acb_mat([[v[i][j,0] for i in mask] for j in range(7)])
        coeff.append(4**len(ids)*gg.det().real*(k+vv*gg.solve(vv.conjugate().transpose())))
    # The positive denominator has no effect on any ordered product's sign.
    def cells(c00,c10,c01,c11,depth):
        faces=[f for f in src["faces"] if set(f)!=set(src["outer"])]
        for f in faces:
            z=acb(1)
            for u,w in zip(f,f[1:]+f[:1]):
                i,j=W.index(u),W.index(w)
                values=[c[i,j] for c in (c00,c10,c01,c11)]
                re=values[0].real
                im=values[0].imag
                for t in values[1:]:
                    re=re.union(t.real)
                    im=im.union(t.imag)
                z*=-acb(re,im)
            if z.imag<0:
                return 1
        assert depth<26,depth
        if depth%2:
            mid0=(c00+c01)/2
            mid1=(c10+c11)/2
            return cells(c00,c10,mid0,mid1,depth+1)+cells(mid0,mid1,c01,c11,depth+1)
        mid0=(c00+c10)/2
        mid1=(c01+c11)/2
        return cells(c00,mid0,c01,mid1,depth+1)+cells(mid0,c10,mid1,c11,depth+1)
    return cells(*coeff,0)


if __name__=="__main__":
    noise=arb(0,arb(".000001"))
    model=fit["model"](p+noise,q+noise)
    _,_,k0=fit["evaluate"](model,[arb(0)]*4,(0,0))
    _,_,ks=fit["evaluate"](model,[1+noise,0,0,0],(0,0))
    assert creation["degrees"](k0)[0]==(0,[0])
    result,gap=creation["degrees"](ks)
    assert result==(3,[3]) and gap>0
    cases=check_pair(0,1)
    print("three unit 724-edge additions create three; cut margin, pair-exclusion cells",gap.str(12),cases)
