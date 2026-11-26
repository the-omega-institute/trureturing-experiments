"""An interval-verified pair of coherent site profiles fitting the 1 TeV quark inputs."""
from pathlib import Path
from runpy import run_path
from flint import acb,acb_mat,arb


m=run_path(str(Path(__file__).with_name("ch39-shell-gw-local-yukawa.py")))
Y=m["local_matrices"]()
adj=m["adjoint"]
eye=m["identity"](3)
up=("87.77892476087915 34.49362739405473", "-101.60703520233089 -41.05288579110566",
    "104.96157181457355 38.03736693474219", "-89.34559663653194 -34.897878611692875",
    "97.27563377404242 40.55466796649334", "-101.07449334237789 -35.76717372137631")
down=("-0.29721454729351665 0.7353702501833489", "0.292741916273106 -0.8672667138658312",
      "-0.22259239433943692 0.8813144906336792", "0.2842694422983706 -0.7584135587507767",
      "-0.3107441583388774 0.833314958131383", "0.22292426273716032 -0.8424774161170383")


def spectrum(weight):
    y=acb_mat(3,3)
    for w,block in zip(weight,Y):
        y+=acb(*[arb(x) for x in w.split()])*block
    h=y*adj(y)
    roots=h.charpoly().roots(tol=1e-32)
    vals=sorted([x.real for x in roots],key=lambda x:float(x.mid()))
    assert all(z.imag.contains(0) and z.real>0 for z in roots)
    proj=[]
    for i,v in enumerate(vals):
        p=eye
        for j,w in enumerate(vals):
            if i!=j:
                p=p*(h-w*eye)/(v-w)
        proj.append(p)
    return [v.sqrt() for v in vals],proj


def trace(a):
    return sum(a[i,i] for i in range(3))


if __name__=="__main__":
    u,p=spectrum(up)
    d,q=spectrum(down)
    targets=[arb(v) for v in (".00000615",".00311",".8616",".0000135",".000268",".01401")]
    mass_error=max(abs(y/t-1).upper() for y,t in zip(u+d,targets))
    vs=[[trace(a*b).real for b in q] for a in p]
    s13=vs[0][2].sqrt();c13=(1-vs[0][2]).sqrt()
    s12=vs[0][1].sqrt()/c13;s23=vs[1][2].sqrt()/c13
    c12=(1-s12*s12).sqrt();c23=(1-s23*s23).sqrt()
    quartet=trace(p[0]*q[1]*p[1]*q[2])
    delta=(quartet+vs[0][2]*s12**2*s23**2*c13**2).arg()
    angles=[x.asin() for x in (s12,s23,s13)]+[delta]
    angle_error=max(abs(a-arb(t)).upper() for a,t in zip(angles,(".22704",".04275",".00377","1.139")))
    assert mass_error<arb("1e-7") and angle_error<arb("1e-7")
    assert arb(".0000315")<quartet.imag<arb(".0000331")
    print("largest relative mass error",mass_error.str(8),"largest absolute angular error",angle_error.str(8))
    print("Jarlskog invariant",quartet.imag.str(14))
    # Complex span contains rank-one matrices by the projective dimension theorem.
    gram=acb_mat([[trace(adj(a)*b) for b in Y] for a in Y])
    assert gram.det().real>0
    total=sum(Y[1:],Y[0])
    for y in Y:
        skew=(adj(total)*y-adj(y)*total)/(2*acb(0,1))
        assert all(m["is_positive"](arb(".002271")*eye+sgn*skew) for sgn in (-1,1))
    alpha=86*arb.pi()/180
    ratio=(arb(".03782")*alpha.cos()-arb(".002271")*alpha.sin())/arb(".34")
    assert ratio>arb(".001")
    print("172-degree profile phase sector mass ratio lower bound",ratio.str(10))
