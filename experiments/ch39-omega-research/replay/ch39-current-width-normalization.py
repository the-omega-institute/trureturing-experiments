"""Recover the absolute normalization from the flavour-inclusive charged-lepton width."""
from pathlib import Path
from runpy import run_path

from flint import arb

root = Path(__file__).parent
inverse = run_path(str(root/"ch39-current-scale-inverse.py"))
four, study, q = inverse["four"], inverse["study"], inverse["q"]


def reference_width(couplings):
    alpha,beta,family = couplings
    charges = study["retuned"]["current"]["weyl_charges"]
    total = arb(0)
    for field, hypercharge in (("L",q("-1/2")), ("ec",q(1))):
        d,r,sign = charges[field]
        assert q(d)/6+q(r)/2 == hypercharge or (q(d)/6+q(r)/2-hypercharge).contains(0)
        total += sum((d*alpha+r*beta+sign*x*family)**2 for x in (1,1,-2))
    return total/(24*arb.pi())


def normalized_lepton_width(rho, eta, c, z):
    unit, colour, right = study["ports"][-3:]
    def a(e, derivative=False):
        return 3*colour(e,unit,derivative)/8 + right(e,unit,derivative)/4

    ceiling = q(".003")
    assert min(c/eta,unit/64) > ceiling
    pole = study["enclose_roots"](lambda e: z+a(e)+4*rho/(c-eta*e), ".0001", ceiling)
    derivative = eta + a(pole,True)*(c-eta*pole)**2/(4*rho)
    assert derivative > 0
    ratio_d = colour(pole,unit)*(c-eta*pole)/(16*rho)
    ratio_r = -right(pole,unit)*(c-eta*pole)/(8*rho)
    multiplicity = 12+27*ratio_d**2+3*(3*ratio_d+ratio_r)**2
    return multiplicity/(24*arb.pi()*derivative)


def recover(shape, width, geometry):
    _, L, rho, eta, c, z = inverse["recover"]([q(1), *shape], geometry)
    f0 = normalized_lepton_width(rho,eta,c,z)/width
    rho, eta, c = rho*f0, eta*f0, c*f0
    target = 15*L*(rho-geometry[-1])
    return target, L, rho, eta, c, z, f0


def verify():
    geometry = inverse["port_moments"]()
    coefficients = inverse["reference_coefficients"](geometry)
    shape = [v/coefficients[0] for v in coefficients[1:]]
    retuned = study["retuned"]
    _, _, (alpha,beta,family), _ = retuned["current"]["light_current"](retuned["logs"])
    width = reference_width((alpha,beta,family))
    time = study["nominal"][2]
    targets = [study["comp"]["model"]["target"](n,time) for n in (11430,11431,11432)]
    exact = recover(shape,width,geometry)
    assert exact[0].overlaps(targets[1]) and abs(exact[0]-targets[1]) < q("1e-12")
    assert exact[-1].overlaps(coefficients[0]) and abs(exact[-1]-coefficients[0]) < q("1e-12")
    shape_error = arb(0,q("1e-9").upper())
    width_error = arb(0,q("1e-5").upper())
    inferred = recover([v*(1+shape_error) for v in shape],width*(1+width_error),geometry)
    target, L, rho, eta, c, z, f0 = inferred
    assert targets[0] > target > targets[2]
    assert abs((f0/coefficients[0]).log()) < q(".0001")
    triplet, octet = four["physical_masses"](c,rho,z)
    thresholds = run_path(str(root/"ch39-neutral-current-threshold-inverse.py"))
    thresholds["exclude_lower_tops"](triplet,octet,c,z/rho,eta,1/rho,time+L.log()/2)
    print("reference inclusive charged-lepton width / mass",width.str(14))
    print("target limits", target.lower().str(15), target.upper().str(15))
    print("absolute log-current normalization error", abs((f0/coefficients[0]).log()).upper().str(10))


if __name__ == "__main__":
    verify()
