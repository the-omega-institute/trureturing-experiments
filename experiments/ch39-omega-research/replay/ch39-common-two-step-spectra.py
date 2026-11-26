"""Four quark ratios on one degree-three Maxwell surface at positive zeta zero 719."""
from pathlib import Path
import runpy
from flint import arb, acb, arb_mat, acb_mat

single=runpy.run_path(str(Path(__file__).with_name("ch39-two-step-created-masses.py")))
fit,flux=single["fit"],single["flux"]
block,wall,J,ii,oo=(fit[k] for k in ("block","W","J","ii","oo"))
rotation=acb(0,arb(".004")).exp()
up=[arb(".00000615")/arb(".8616"),arb(".00311")/arb(".8616")]
down=single["target"]
height=acb.zeta_zero(719).imag
angles=[fit["s"][key]*height for key in ("a","b")]


def constraints(model,weights,goals,columns):
    m,orbits=model
    m=acb_mat(m)
    for val,orbit in zip(weights,orbits):
        for h in orbit:
            m+=val*h*h.conjugate().transpose()
    proj=block(m,wall,J)*block(m,J,J).inv()
    k=block(m,wall,wall)-proj*block(m,J,wall)
    eig=fit["symbols"](-block(k,ii,oo))
    def top(mat):
        return -mat[wall.index(7),wall.index(6)]*rotation
    residual=[abs(eig[j+1])**2-goals[j]*abs(eig[0])**2 for j in range(2)]+[top(k).imag]
    columns_out=[]
    for j in columns:
        deriv=acb_mat(7,7)
        for h in orbits[j]:
            v=block(h,wall,(0,))-proj*block(h,J,(0,))
            deriv+=v*v.conjugate().transpose()
        de=fit["symbols"](-block(deriv,ii,oo))
        columns_out.append([2*(eig[l+1].conjugate()*de[l+1]-goals[l]*eig[0].conjugate()*de[0]).real
                            for l in range(2)]+[top(deriv).imag])
    return arb_mat([[v] for v in residual]),arb_mat(columns_out).transpose(),k,eig,top(k)


def solve(model,center,goals,columns):
    center=list(map(arb,center))
    for _ in range(7):
        r,j,*_=constraints(model,center,goals,columns)
        delta=j.solve(r)
        for t,i in enumerate(columns):
            center[i]=(center[i]-delta[t,0]).mid()
    # The centers are Newton trial points. The solution lies in this entire box.
    radius=arb("1e-18")
    box=[v+arb(0,radius) if i in columns else v for i,v in enumerate(center)]
    angular=fit["model"](*(v+arb(0,arb("1e-24")) for v in angles))
    inv=constraints(model,center,goals,columns)[1].inv()
    r=constraints(angular,center,goals,columns)[0]
    _,j,k,e,x=constraints(angular,box,goals,columns)
    step=inv*r
    err=arb_mat([[int(i==j) for j in range(3)] for i in range(3)])-inv*j
    for i in range(3):
        assert abs(step[i,0])+radius*sum(abs(err[i,v]) for v in range(3))<radius
    assert min(box)>0 and x.real<0 and abs(e[0])>1
    assert flux["degrees"](k)[0]==(3,[3])
    print("four clock strengths",*(v.str(12) for v in box))
    return k


def shared_surface():
    model=fit["model"](*angles)
    ku=solve(model,("3.690","3.60206356","2.70565405","2.73355298"),up,(1,2,3))
    kd=solve(model,("1.57063427",".877253274",".0084002227",".00004"),down,(0,1,2))
    i,j=wall.index(11),wall.index(12)
    yu,yd=-ku[i,j],-kd[i,j]
    bare=next(z for a,b,d,z in fit["s"]["edges"] if (a,b,d)==(11,12,392))
    # Im[(Y_d+w z) conj(Y_u)]=0. Positivity fixes the same ray.
    w=-(yd*yu.conjugate()).imag/(bare*yu.conjugate()).imag
    assert 0<w<arb(".04") and ((yd+w*bare)*yu.conjugate()).real>0
    for a,b in ((11,12),(12,18),(18,11)):
        bond=bare
        at,bt=wall.index(a),wall.index(b)
        kd[at,at]+=w
        kd[bt,bt]+=w
        kd[at,bt]-=w*bond
        kd[bt,at]-=w*bond.conjugate()
    for k in (ku,kd):
        observed,margin=flux["degrees"](k)
        assert observed==(3,[3]) and margin>arb(".0001")
    bare_k=constraints(model,[arb(0)]*4,down,(0,1,2))[2]
    assert flux["degrees"](bare_k)[0]==(0,[0])
    print("added down-type stiffness on the 392 bottom orbit",w.str(12))
    return ku,kd


if __name__=="__main__":
    shared_surface()
