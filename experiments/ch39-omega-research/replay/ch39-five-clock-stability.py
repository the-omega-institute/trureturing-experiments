"""Local response, near-minimum arithmetic density, and stability of inverse design."""
from pathlib import Path
import runpy
from flint import arb,arb_mat

p=runpy.run_path(str(Path(__file__).with_name("ch39-five-clock-selection.py")))
q,pi=p["q"],p["pi"]
K=q("9e-8")
G0=q("1.6150e-10")
G2=q("3.7360e-8")
target=p["target"]
w=p["weights"]
C=p["C"]
H0=p["eval_cost"](*target,2)[3:]
D=q("1.0351e-10")+q("2.2182e-11")
pert=G2+80*D
def extremal(a,b,c):
    root=((a-c)**2+4*b*b).sqrt()
    return (a+c-root)/2,(a+c+root)/2

def check():
    lo,hi=extremal(*H0)
    a,b=K*lo,K*hi
    assert min(a.lower(),b.lower())>pert
    rmin=(a-pert)*(b-pert)
    rmax=(a+pert)*(b+pert)
    slope=(1/(pi*rmax.sqrt())).lower().union((1/(pi*rmin.sqrt())).upper())
    shift=q("1e-6")
    pts=[x+arb(0,shift.upper()) for x in target]
    rows=p["data"](*pts,2)
    Cx=arb_mat([[rows[j][i] for j in (1,2)] for i in (1,2)])
    Ci=Cx.inv()
    norm=abs(Ci[0,0]).upper()+abs(Ci[1,0]).upper()+abs(Ci[0,1]).upper()+abs(Ci[1,1]).upper()
    hs=[sum(z*k for z,k in zip(col,w)) for col in zip(*rows)][3:]
    H=max((abs(hs[0])+abs(hs[1])).upper(),(abs(hs[1])+abs(hs[2])).upper())
    # Modify the two weights from the fully compensated 719 value. Its gradient is zero.
    # The Euclidean spectral Hessian cap converts to row-sum norm by <=sqrt(2).
    gdiff=(K*H+q(2).sqrt()*pert)*shift
    extra=(norm*gdiff).upper()
    ext=K*q(".015")-G0-9*(D+extra)
    loc=5*K-G2-80*(D+extra)
    assert Cx.det()>0 and ext>0 and loc>0
    assert min(w[i].lower() for i in (1,2))*K>D+extra
    print("clock Hessian eigenvalue enclosures",lo.str(12),hi.str(12))
    print("density slope",slope.lower().str(12),slope.upper().str(12))
    print("target square radius",shift.str(12),"extra weight bound",extra.str(12),
          "exterior margin",ext.lower().str(10),"local curvature margin",loc.lower().str(10))
if __name__=="__main__":
    check()
