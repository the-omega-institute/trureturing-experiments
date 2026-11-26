"""Check the quark-mass determinant's exact dependence on the CKM phase."""
from pathlib import Path
from runpy import run_path
from flint import acb, acb_mat, arb

m=run_path(str(Path(__file__).with_name("ch39-shell-gw-195-flavor.py")))
eye,adj=m["eye"],m["adj"]
up=list(map(arb,(".00000615",".00311",".8616")))
down=list(map(arb,(".0000135",".000268",".01401")))
angles=list(map(arb,(".22704",".00377",".04275")))


def ckm(delta):
    result=eye(3)
    for i,j,t,phase in ((1,2,angles[2],0),(0,2,angles[1],-delta),(0,1,angles[0],0)):
        U=eye(3)
        U[i,i]=U[j,j]=t.cos()
        U[i,j]=t.sin()*acb(0,phase).exp()
        U[j,i]=-U[i,j].conjugate()
        result=result*U
    return result


def determinant(delta,shift,a,b):
    U=acb_mat([[a*up[i]**2 if i==j else 0 for j in range(3)] for i in range(3)])
    D=acb_mat([[b*down[i]**2 if i==j else 0 for j in range(3)] for i in range(3)])
    V=ckm(delta)
    return (shift*eye(3)+U+V*D*adj(V)).det().real


if __name__=="__main__":
    t12,t13,t23=angles
    k=2*t12.sin()*t12.cos()*t23.sin()*t23.cos()*t13.sin()
    beta=k*(up[2]**2-up[1]**2)*(down[1]**2-down[0]**2)
    assert beta>0
    measured=arb("1.139")
    for shift,a,b in ((1,1,1),(arb(".01"),2,7),(0,13,3)):
        factor=beta*a*b*(shift+a*up[0]**2+b*down[2]**2)
        diff=determinant(measured,shift,a,b)-determinant(0,shift,a,b)
        assert (diff-factor*(1-measured.cos())).contains(0) and diff>0
    diff=determinant(measured,1,1,1).log()-determinant(0,1,1,1).log()
    print("phase coefficient",beta.str(12),"central logarithmic cost above CP-even minimum",diff.str(12))
