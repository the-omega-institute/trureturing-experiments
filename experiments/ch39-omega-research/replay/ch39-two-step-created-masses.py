"""A positive degree-three two-step down-type spectrum at fixed zeta heights."""
from pathlib import Path
import runpy
from flint import arb, acb, acb_mat, arb_mat

flux=runpy.run_path(str(Path(__file__).with_name("ch39-source-created-flux.py")))
fit=flux["fit"]
block,wall,J,ii,oo=(fit[k] for k in ("block","W","J","ii","oo"))
target=[arb(".000268")/arb(".01401"),arb(".0000135")/arb(".01401")]


def system(model,weights,columns):
    m,orbits=model
    m=acb_mat(m)
    for val,orbit in zip(weights,orbits):
        for h in orbit:
            m+=val*h*h.conjugate().transpose()
    proj=block(m,wall,J)*block(m,J,J).inv()
    k=block(m,wall,wall)-proj*block(m,J,wall)
    symbols=fit["symbols"](-block(k,ii,oo))
    eq=[abs(e)**2-r*abs(symbols[0])**2 for e,r in zip(symbols[1:],target)]
    deriv=[]
    for j in columns:
        t=acb_mat(3,3)
        for h in orbits[j]:
            v=block(h,wall,(0,))-proj*block(h,J,(0,))
            t-=block(v,ii,(0,))*block(v,oo,(0,)).conjugate().transpose()
        ds=fit["symbols"](t)
        deriv.append([2*(symbols[i+1].conjugate()*ds[i+1]-r*symbols[0].conjugate()*ds[0]).real
                      for i,r in enumerate(target)])
    return arb_mat([[v] for v in eq]),arb_mat(deriv).transpose(),k,symbols


def certify(zero,columns,start):
    height=acb.zeta_zero(zero).imag
    angles=[fit["s"][key]*height for key in ("a","b")]
    exact=fit["model"](*angles)
    center=[arb(0)]*4
    for i,v in zip(columns,start):
        center[i]=arb(v)
    for _ in range(6):
        res,j,_,_=system(exact,center,columns)
        shift=j.solve(res)
        for i,col in enumerate(columns):
            center[col]=(center[col]-shift[i,0]).mid()
    radius=arb("1e-6")
    angular=fit["model"](*(a+arb(0,arb("1e-10")) for a in angles))
    box=[v+arb(0,radius) if j in columns else v for j,v in enumerate(center)]
    assert all(box[j]>0 for j in columns)
    r,_,_,_=system(angular,center,columns)
    inv=system(exact,center,columns)[1].inv()
    _,Jbox,k,e=system(angular,box,columns)
    err=arb_mat([[1,0],[0,1]])-inv*Jbox
    step=inv*r
    for i in range(2):
        assert abs(step[i,0])+radius*sum(abs(err[i,j]) for j in range(2))<radius
    assert 0<abs(e[2])<abs(e[1])<abs(e[0])
    degree,margin=flux["degrees"](k)
    print(zero,*(v.str(12) for v in center),degree,margin.str(10))
    assert degree==(3,[3]) and margin>arb(".0001")
    assert flux["degrees"](system(angular,[arb(0)]*4,columns)[2])[0]==(0,[0])


def rounded_neighborhood(zero,weights,columns):
    height=acb.zeta_zero(zero).imag
    noise=arb(0,arb(".0001"))
    angular=fit["model"](*(fit["s"][key]*height+noise for key in ("a","b")))
    weights=[arb(w)+noise if i in columns else arb(0) for i,w in enumerate(weights)]
    _,_,k,e=system(angular,weights,columns)
    assert flux["degrees"](k)[0]==(3,[3])
    assert flux["degrees"](system(angular,[arb(0)]*4,columns)[2])[0]==(0,[0])
    assert arb(".000265")/arb(".01410")<abs(e[1]/e[0])**2<arb(".000271")/arb(".01392")
    assert arb(".0000133")/arb(".01410")<abs(e[2]/e[0])**2<arb(".0000137")/arb(".01392")


if __name__=="__main__":
    certify(8,(2,3),("1.48","1.29"))
    certify(285,(0,1),("1.57",".85"))
    rounded_neighborhood(8,(0,0,"1.477","1.293"),(2,3))
    rounded_neighborhood(285,("1.573",".854",0,0),(0,1))
