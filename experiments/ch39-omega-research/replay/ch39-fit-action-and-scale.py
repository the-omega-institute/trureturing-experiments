"""Dynamical cost and the intersector mass normalization of the complete finite CKM fit."""
from pathlib import Path
import runpy
from flint import arb,acb_mat

metric=runpy.run_path(str(Path(__file__).with_name("ch39-full-metric-ckm.py")))
fit,finite=metric["fit"],metric["finite"]
adj=metric["adj"]
bonds=acb_mat(19,6)
for j in range(6):
    i,k=fit["inner"][j%3],fit["outer"][(j%3+j//3+1)%3]
    bonds[i,j]=1
    bonds[k,j]=-finite["bond"](i,k).conjugate()


def sector(q,weights,up):
    m=metric["model"][0]+sum((w*v for w,v in zip(weights[:4],metric["source"])),acb_mat(19,19))
    W,J=fit["W"],fit["J"]
    block=fit["block"]
    k=block(m,W,W)-block(m,W,J)*block(m,J,J).solve(block(m,J,W))
    for cycle,angle in (((7,6,8),arb.pi()-arb(".002")),((11,12,18),arb(".03606860"))):
        w=finite["fix_phase"](k,cycle,angle)
        for i,j in zip(cycle,cycle[1:]+cycle[:1]):
            z=finite["bond"](i,j)
            m[i,i]+=w
            m[j,j]+=w
            m[i,j]-=w*z
            m[j,i]-=w*z.conjugate()
    r=adj(bonds)*m.solve(bonds)
    assert min(abs(r[0,3]).lower()-abs(r[0,j]).upper() for j in (4,5))>arb(".075")
    w=weights[4:]
    totals=[sum(w[:3]),sum(w[3:])]
    pure=[totals[0],arb(0),arb(0),totals[1],arb(0),arb(0)]
    def action(v):
        return acb_mat([[int(i==j)+v[i]*r[i,j] for j in range(6)] for i in range(6)]).det().real
    def surface(v):
        mat=acb_mat(k)
        for j,val in enumerate(v):
            finite["add_edge"](mat,fit["inner"][j%3],fit["outer"][(j%3+j//3+1)%3],val)
        return mat
    if not up:
        response=r*acb_mat([[int(i==j)+w[i]*r[i,j] for j in range(6)] for i in range(6)]).inv()
        # Move load from the larger second negative-orientation arc to the smaller first one.
        assert w[4]-w[3]>arb(".28")
        assert response[4,4].real-response[3,3].real>arb(".0135")
        exchange=list(w)
        exchange[3]+=arb(".01")
        exchange[4]-=arb(".01")
        assert finite["flux"]["degrees"](surface(exchange))[0]==(3,[3])
        assert action(w)/action(exchange)>arb("1.000135")
    ratio=action(w)/action(pure)
    assert ratio>(arb("1.019") if up else arb("1.026"))
    for j in range(3):
        other=[totals[0],arb(0),arb(0),arb(0),arb(0),arb(0)]
        other[3+j]=totals[1]
        assert finite["flux"]["degrees"](surface(other))[0]==((3,[1,2,3]) if up else (3,[3]))
    nearby=[(9*a+b)/10 for a,b in zip(w,pure)]
    assert finite["flux"]["degrees"](surface(nearby))[0]==(3,[3])
    assert action(w)/action(nearby)>(arb("1.0011") if up else arb("1.0013"))
    t=-block(surface(w),fit["ii"],fit["oo"])
    D=acb_mat([[q,0,0],[0,1,0],[0,0,1]])
    y=D*adj(metric["F"])*(t*t)*metric["F"]*D
    return max(z.real for z in (y*adj(y)).eig(algorithm="rump")).sqrt()


if __name__=="__main__":
    dilution=(1+metric["step"]/(2*metric["qlo"]))**2
    for i in range(metric["ncell"]):
        q=metric["qlo"]+(arb(i)+arb(1)/2)*metric["step"]
        w=metric["weight_row"](q)
        su,sd=sector(q,w[:10],True),sector(q,w[10:],False)
        ratio=sd/su
        assert arb("1.082")<ratio/dilution<ratio*dilution<arb("1.093")
        goal=arb(".01401")/arb(".8616")
        assert arb(".01487")<goal/(ratio*dilution)<dilution*goal/ratio<arb(".01504")
    print("all metric cells: top/bottom scales require independent relative factor in (.01487,.01504)")
