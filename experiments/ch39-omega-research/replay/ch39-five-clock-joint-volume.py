"""Pay for the selector and the physical neutral vector sector in the same flat-space energy equation."""
from fractions import Fraction as F
from pathlib import Path
import runpy
from flint import arb

root=Path(__file__).resolve().parent
pent=runpy.run_path(str(root/"ch39-five-clock-selection.py"))
cross=runpy.run_path(str(root/"ch39-spin-gauge-crossing.py"))
short=pent["m"]["short"]
q=pent["q"]
K=q("9e-8")
# Fixed target, arbitrary geometry/source in the proven compact box. The physical
# neutral-field realization leaves exactly three adjoint polarization differences.
A=pent["station"]
coef=[q("1.0351e-10"),q("2.2182e-11")]
portal=sum(x*A[i][0] for i,x in zip((1,2),coef)).upper()
cl=K*pent["emin"]+arb(0,portal)
neutral=3*q("8.37729e-12").union(q("9.37987e-12"))
add=cl+neutral
desired=-2*arb.pi()**2*add
assert desired<0 and abs(desired)<q(".024")

def matter():
    # T = triplet heat; S = scalar graph heat. The vector content is T+S.
    # For a massless vector one extra minimally coupled real scalar restores
    # its massive Proca determinant in the zero-mass limit, with gauge zero modes omitted.
    ns,nv,nw=20,4,16
    assert ns+3*nv-2*nw==0
    assert ns*(F(1,6)-F(1,5))-nv*F(1,2)+nw*F(1,6)==0
    assert F(ns,180)-F(nv,15)+F(7*nw,720)==0
    assert -F(ns,180)+F(29*nv,60)+F(nw,90)==2

def faces():
    # The three projected spectral residuals may use this continuous extra
    # vacuum energy. The checked volume faces allow all targets up to .024.
    g=cross["b"]["model"]["gauge_inputs"]
    for lh,rh,hmax in ((144,4716,30),(216,4068,27)):
        _,z,_,_=cross["fit"](lh,rh,hmax,g)
        shift=cross["b"]["shift"]
        assert 18*(short["s0"]+shift.upper())<z<18*(short["s1"]+shift.lower())

def curvature():
    gr=runpy.run_path(str(root/"ch39-short-arc-curvature.py"))
    err=gr["error"]()
    s0,s1=short["s0"],short["s1"]
    stack=[([short["full"]]*3,s0.union(s1),0)]
    vals=[]
    while stack:
        box,s,depth=stack.pop()
        val=gr["assess"](box,s,err)
        if val>q(".05") and val<2:
            vals.append(val)
            continue
        if any(not (short["margin"](i,1,box,s)-short["P"][i][2]*desired).contains(0) for i in range(3)):
            continue
        assert depth<24
        j=depth%4
        ax=s if j==3 else box[j]
        mid=ax.mid()
        for sub in (ax.lower().union(mid),mid.union(ax.upper())):
            child=box[:]
            if j<3:
                child[j]=sub
            stack.append((child,sub if j==3 else s,depth+1))
    assert vals
    print("charged curvature moment",min(v.lower() for v in vals).str(12),max(v.upper() for v in vals).str(12))

if __name__=="__main__":
    matter()
    faces()
    curvature()
    print("classical selector energy / Q^4",cl.lower().str(12),cl.upper().str(12))
    print("required charged I2",desired.lower().str(12),desired.upper().str(12))
