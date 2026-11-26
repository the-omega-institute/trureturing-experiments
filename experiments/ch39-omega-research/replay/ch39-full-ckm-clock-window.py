"""An open rectangle of two RH angles admits the common-scale finite fit."""
from pathlib import Path
import runpy
from flint import arb,acb_mat

scaled=runpy.run_path(str(Path(__file__).with_name("ch39-common-yukawa-scale.py")))
metric,fit,finite=(scaled[k] for k in ("metric","fit","finite"))
radius=arb("2e-7")
angles=metric["common"]["angles"]
noise=arb(0,radius)
vary=fit["model"](*(v+noise for v in angles))
base=metric["model"]


def square(k,weights):
    for i in range(3):
        finite["add_edge"](k,fit["inner"][i],fit["outer"][i],weights[10])
        for end in (1,2):
            finite["add_edge"](k,fit["inner"][i],fit["outer"][(i+end)%3],weights[4+3*(end-1)+i])
    finite["fix_phase"](k,(7,6,8),arb.pi()-arb(".002"))
    finite["fix_phase"](k,(11,12,18),arb(".03606860"))
    degree,margin=finite["flux"]["degrees"](k)
    assert degree==(3,[3]) and margin>arb(".0001")
    return -fit["block"](k,fit["ii"],fit["oo"])


def generators(model):
    m,orbits=model
    return m,[sum((v*metric["adj"](v) for v in orb),acb_mat(19,19)) for orb in orbits]


def transfer(model,w):
    m,gen=model
    m=m+sum((weight*col for weight,col in zip(w,gen)),acb_mat(19,19))
    block,W,J=fit["block"],fit["W"],fit["J"]
    y=block(m,W,W)-block(m,W,J)*block(m,J,J).solve(block(m,J,W))
    return square(y,w)


def frobenius(m):
    return sum(abs(m[i,j]).upper()**2 for i in range(m.nrows()) for j in range(m.ncols())).sqrt()


def perturbation(center,varying,q,w):
    t=transfer(center,w)
    change=transfer(varying,w)-t
    f=metric["F"]
    D=acb_mat([[q,0,0],[0,1,0],[0,0,1]])
    y=D*metric["adj"](f)*(t*t)*f*D
    # D is a contraction; changing a ring matrix changes its square in three terms.
    bound=2*frobenius(t)*frobenius(change)+frobenius(change)**2
    heavy=max(v.real for v in (y*metric["adj"](y)).eig(algorithm="rump")).sqrt()
    return bound/heavy


def conclude(max_err):
    # The metric dilation also bounds the relative matrix error at the actual q.
    upper=[v*(1+metric["step"]/(2*metric["qlo"]))**2 for v in max_err]
    assert upper[0]<arb("2.4e-8") and upper[1]<arb("1.8e-7")
    gaps=[arb(".0035"),arb(".018")]
    dv=sum(2*err/(gap-err) for err,gap in zip(upper,gaps))
    targets=[metric["common"]["up"],list(reversed(metric["common"]["down"]))]
    for err,rs in zip(upper,targets):
        for ref in rs:
            a=ref*(1+arb(0,arb(".00013")))+arb(0,err)
            b=1+arb(0,err)
            assert abs((a/b)/ref-1)<arb(".0036")
    target=arb(".01401")/arb(".8616")
    ratio=target*(1+arb(0,arb(".00013")))*(1+arb(0,upper[1]))/(1+arb(0,upper[0]))
    assert abs(ratio/target-1)<arb(".000131")
    # Compare to the physical-q reference, whose angular and phase errors were certified earlier.
    a,b,c=[arb(x)+arb(0,arb(".000019")) for x in (".22501",".04183",".003732")]
    cc=(1-c*c).sqrt()
    nc=c+arb(0,dv)
    na=(a*cc+arb(0,dv))/(1-nc*nc).sqrt()
    nb=(b*cc+arb(0,dv))/(1-nc*nc).sqrt()
    assert all(abs(v-arb(t))<arb(".000055") for v,t in zip((na,nb,nc),(".22501",".04183",".003732")))
    abscenter=[(a*cc).lower(),(b*cc).lower(),c.lower(),
               (((1-a*a)*(1-b*b)).sqrt()-a*b*c).lower()]
    def offset(l,m,n):
        return (l*m*n/((1-l*l)*(1-m*m)).sqrt()).asin()
    # Endpoint phase-difference bound from disks about the four reference elements.
    extra=sum((dv/x).asin() for x in abscenter)+offset(a,b,c)+offset(na,nb,nc)
    delta=arb("1.147")+arb(0,(arb(".006")+extra).upper())
    assert abs(delta-arb("1.147"))<arb(".017")
    j=na*nb*nc*((1-na*na)*(1-nb*nb)).sqrt()*(1-nc*nc)*delta.sin()
    assert arb(".0000300")<j<arb(".0000325")
    density=(radius/arb.pi())**2
    assert density>arb("4.05e-15")
    print("clock square density, extra CKM entry error, J",density.str(10),dv.upper().str(9),j.str(8))


if __name__=="__main__":
    central,varying=generators(base),generators(vary)
    max_err=[arb(0),arb(0)]
    for i in range(metric["ncell"]):
        q=metric["qlo"]+(arb(i)+arb(1)/2)*metric["step"]
        w=scaled["weights"](q)
        for j,part in enumerate((w[:11],w[11:])):
            e=perturbation(central,varying,q,part).upper()
            assert e.is_finite() and e>0
            max_err[j]=max(max_err[j],e)
    print("relative angular matrix error upper bounds",*(v.str(10) for v in max_err))
    conclude(max_err)
