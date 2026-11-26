"""All-direction initial adjoint screening with a fundamental and an exterior-square Dirac field."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

ROOT = Path(__file__).parent
MATCH = run_path(str(ROOT/"ch39-quartic-soft-adjoint-match.py"))
a = MATCH["a"]


def basin():
    # Rescaled RG time, with z=g5^2 and xi=|y|^2/z.
    step = a("1/4096")
    margin = 16*(-step).union(step)
    five,ten = a("7/2"),a("13/10")
    l5,l10 = a("206/15"),a("314/15")
    def rhs(f,t):
        return f*(a("16/5")*f+3*t-l5),t*(f+a("44/5")*t-l10)
    for j in range(2048):
        if a("16/5")*five.upper()+3*ten.upper()<l5 and five.upper()+a("44/5")*ten.upper()<l10:
            print("Invariant decreasing rectangle by scaled time",(j*step).str(7))
            return
        df,dt = rhs(five+margin,ten+margin)
        assert abs(df)<16 and abs(dt)<16
        five,ten = five+step*df,ten+step*dt
    raise AssertionError("Downward invariant rectangle not reached")


def small_representations():
    rep = run_path(str(ROOT/"ch39-quartic-adjoint-lock.py"))
    for labels,dim,_ in rep["small"](5):
        if dim==1:
            continue
        casimir = sum(F(min(i,j)*(5-max(i,j)),10)*labels[i-1]*(labels[j-1]+2)
                      for i,j in product(range(1,5),repeat=2))
        assert F(dim,24)*casimir=={5:F(1,2),10:F(3,2),15:F(7,2)}[dim]
    # For dimension >=24, a non-fundamental has Casimir at least 18/5.
    assert F(24,24)*F(18,5)>F(5,2)


def balanced_masses(r, z):
    y5,y10 = (a("7/2")*z).sqrt(),(a("13/10")*z).sqrt()
    d = y10*r
    real = ((1+d*d/2)/2).sqrt()
    imag = (1-real**2).sqrt()
    c = ((real+d)**2+imag**2).sqrt()
    e = ((real-d)**2+imag**2).sqrt()
    assert (c*e-1).contains(0)
    gap = y5*r
    # The weak fundamental mass has modulus 1. The color mass has modulus 1/c.
    wreal = (1/c**2-1-gap**2)/(2*gap)
    assert abs(wreal)<1
    wimag = (1-wreal**2).sqrt()
    color = ((wreal+gap)**2+wimag**2).sqrt()
    assert (color*c-1).contains(0)
    # C,1,E are the three masses of the exterior square, all in units of v.
    masses = [color,arb(1),c,arb(1),e]
    indices = [(F(1,5),0,F(1,2)),(F(3,10),F(1,2),0),
               (F(4,5),0,F(1,2)),(F(1,10),F(3,2),1),(F(3,5),0,0)]
    for i in range(3):
        assert sum(a(idx[i])*m.log() for idx,m in zip(indices,masses)).contains(0)
    return c


def verify():
    ctx.prec = 140
    small_representations()
    basin()
    e, f = a("1/3072"),a("1/12")
    xi5,xi10 = a("7/2"),a("13/10")
    for label,kind,r,v,state in MATCH["matched_points"]():
        z = 1/state[0]
        # Tr_10(L(Phi)^4) = -3 Tr Phi^4 + 3 (Tr Phi^2)^2.
        a4 = (15-2*xi5**2+6*xi10**2)*z*z
        a4 += (-60+4*(xi5+3*xi10))*z*e+a("64/5")*e**2-30*f**2
        a22 = (9-6*xi10**2)*z*z+a("168/25")*e**2+30*f**2
        assert a22 < 0 and a4+a22 < 0
        # A single exterior-square Dirac cannot screen diag(4,-1,-1,-1,-1).
        shape = a("13/20")
        limit = a("76/33")
        lower = (15*shape+9-6*(1-shape)*limit**2)*z*z-60*shape*z*e
        assert lower > 0
        c = balanced_masses(r,z)
        print(label,kind,"upper endpoints",a22.str(7),(a4+a22).str(7),"m_cc/v",c.str(8))


if __name__ == "__main__":
    verify()
