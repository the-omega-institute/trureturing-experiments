"""Sharp mass-phase bounds when all three real one-loop thresholds cancel."""
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

ROOT = Path(__file__).parent
MIX = run_path(str(ROOT/"ch39-mixed-adjoint-screening.py"))
a, MATCH = MIX["a"],MIX["MATCH"]


def phase_bounds(r,z):
    d = (a("13/10")*z).sqrt()*r
    h = (a("7/2")*z).sqrt()*r
    def stretch(t):
        v = d/t
        return (1+v*v+v*(2+v*v).sqrt()).sqrt()
    def endpoint(fn,lo,hi):
        left,right = F(lo),F(hi)
        assert fn(a(left)) < 0 < fn(a(right))
        for _ in range(36):
            mid = (left+right)/2
            value = fn(a(mid))
            if value < 0:
                left = mid
            elif value > 0:
                right = mid
            else:
                break
        assert right-left < F(1,10**8)
        return a(left).union(a(right))
    tmin = endpoint(lambda t: h*t**3-1+1/stretch(t),".4","1.1")
    tmax = endpoint(lambda t: h*t**3-1-stretch(t),"1.1","2.5")
    def center(t):
        real = ((t*t+d*d/2)/2).sqrt()
        return acb(real,(t*t-real*real).sqrt())
    qmin = center(tmin)
    sine = qmin.imag/abs(qmin+d/5)
    chord = (d/tmax)**2
    assert a(".6") < sine < a(".8") and a(".008") < chord < a(".02")

    def height(t):
        return ((t/d)**4-a("1/4")).sqrt()
    assert abs(1-1/stretch(arb(1))) < h < 1+1/stretch(arb(1))
    assert abs(1-stretch(arb(1))) < h < 1+stretch(arb(1))
    bounds = [height(tmin),height(tmax)]
    assert [b.zeta_nzeros().unique_fmpz() for b in bounds] == [1,23]
    # Endpoints are strictly between zeros by interval counting on expanded boxes.
    for b,count in zip(bounds,(1,23)):
        assert (b+a(".001")*(-arb(1)).union(arb(1))).zeta_nzeros().unique_fmpz()==count
    negative = endpoint(lambda t: h*t**4/stretch(t)-d,".6","1.2")
    # A common nonabelian determinant phase is attainable on the positive-real-part branch.
    t = endpoint(lambda t: h*t**4*stretch(t)-d,".6","1.2")
    D = center(t)
    C,E = D+d,D-d
    W,P = -(h/d)*C,-(h/d)*D
    assert (P-W-h).contains(0)
    assert (abs(W*D**3)-1).contains(0)
    assert (P*C*D**2-W*D**3).contains(0)
    assert (abs(C*E/D**2)-1).contains(0)
    assert abs(C*E/D**2-1)>a(".008")
    for locked in (t,negative):
        assert a("17.8") < height(locked) < a("21")
    print("RH range",bounds[0].str(8),bounds[1].str(8),
          "aligned heights",height(t).str(8),height(negative).str(8))
    return sine,chord,t


def verify():
    ctx.prec = 150
    assert acb.zeta_zero(1).imag < a('17.8') and acb.zeta_zero(2).imag > a('21')
    for n,kind,r,v,state in MATCH["matched_points"]():
        sine,chord,t = phase_bounds(r,1/state[0])
        print(n,kind,"minimum sine",sine.str(9),"minimum phase chord",chord.str(9),
              "aligned nonabelian phases at |D|",t.str(9))


if __name__ == "__main__":
    verify()
