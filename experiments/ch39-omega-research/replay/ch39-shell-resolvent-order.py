"""Rational all-scale ordering and an explicit heat-trace crossing for two Berry shells."""
from pathlib import Path
from runpy import run_path
from flint import arb, fmpq, fmpq_poly


local = run_path(str(Path(__file__).with_name("ch39-shell-berry-rh-matching.py")))
spectra = {r:local["eigenvalue_boxes"](local["covariant_laplacian"](r,3)) for r in (17,29)}
X = fmpq_poly([0,1])


def characteristic(r):
    result = fmpq_poly([1])
    for lo,hi,n in spectra[r]:
        cut = lo if r == 17 else hi
        numerator = (cut*10**12).unique_fmpz()
        assert numerator is not None
        shift = fmpq(12,10**6) if r == 29 else 0
        result *= (X + fmpq(numerator,10**12) + shift)**n
    return result


def heat(r,y):
    return sum(n*(-y*lo.union(hi)).exp() for lo,hi,n in spectra[r])


if __name__ == "__main__":
    a,b = characteristic(17),characteristic(29)
    cross = a*b.derivative() - b*a.derivative()
    coeff = cross.coeffs()
    assert len(coeff)==21 and all(c>0 for c in coeff)
    print("strictly positive exact rational coefficients of the resolvent difference",len(coeff))
    rhs = local["target"]
    trace_lower = sum(n*local["Z"]((17+10**6*hi).sqrt()) for _,hi,n in spectra[17])/local["Z"](arb(0))
    assert trace_lower>rhs
    print("17-shell zero-scale trace at stiffness one million exceeds",trace_lower.lower().str(12))

    def heat_gap(y):
        x = (heat(17,y)/rhs).log()/17
        assert x>0
        return (-29*x).exp()*heat(29,y)-rhs
    lo,hi = arb(182),arb(185)
    assert heat_gap(lo)<0 and heat_gap(hi)>0
    for _ in range(20):
        mid = (lo+hi)/2
        val = heat_gap(mid)
        if val<0:
            lo = mid
        elif val>0:
            hi = mid
        else:
            break
    assert heat_gap(lo)<0 and heat_gap(hi)>0
    y = lo.union(hi)
    x = (heat(17,y)/rhs).log()/17
    print("common heat matching y",y.str(9),"Q/mu",x.rsqrt().str(9),"stiffness",(y/x).str(9))
