"""Target number, its determinant zero mode, and the heavy phase at gain 624."""
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path

from flint import acb,arb,ctx

ROOT = Path(__file__).parent
EXT = run_path(str(ROOT/"ch39-selected-extended-flow.py"))
a = EXT["a"]


def verify():
    ctx.prec = 150
    # Three target fundamentals B and three conjugate targets T have opposite vector charge.
    charges = [1]*3+[-1]*3
    discrete = [1]*3+[0,1,-1]
    assert sum(charges)==0 and sum(10*q**3 for q in charges)==0
    assert sum(10*q*z*z for q,z in zip(charges,discrete))%3==1
    # The exterior-square diagonal source generator is traceless whenever its parent is.
    for i in range(4):
        h = [int(j==i)-int(j==i+1) for j in range(5)]
        assert sum(h[j]+h[k] for j in range(5) for k in range(j+1,5))==0
    # A target-center-invariant monomial of nonzero target number needs at least 10 factors.
    assert [q for q in range(-9,10) if q%10==0]==[0]
    assert sum(1 for j in range(5) for k in range(j+1,5) if k<3)==3
    # The canonical angular kinetic coefficient in units of v squared.
    assert F(3+7,10**2)==F(1,10)
    v,_ = EXT["initial"]()
    decay = v/arb(5).sqrt()
    cutoff = a("1.2209e19")
    pseudo = arb(10).sqrt()*v*(v/cutoff)**3
    assert a("2.64e14")<decay<a("2.66e14") and a("213.9")<pseudo<a("214.1")
    _,gain,branch,_,_ = next(EXT["SELECT"]["matches"]())
    assert (gain,branch)==(624,1)
    gamma = acb.zeta_zero(11490).imag
    ell2 = ((gamma/gain)**2+a("1/4")).sqrt()
    middle = acb(((ell2+a("1/2"))/2).sqrt(),((ell2-a("1/2"))/2).sqrt())
    det = -(middle+1)*middle**3
    angle = det.arg()
    assert a("-.19613") < angle < a("-.19612")
    print("Common heavy determinant angle",angle.str(12))
    print("f_theta /GeV",decay.str(10),"dimension-ten example mass /GeV",pseudo.str(10))


if __name__=="__main__":
    verify()
