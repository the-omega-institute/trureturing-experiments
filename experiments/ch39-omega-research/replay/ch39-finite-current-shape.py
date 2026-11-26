"""Enclose four low-energy shape coefficients from finite noisy Euclidean samples."""
from math import prod
from pathlib import Path
from runpy import run_path

from flint import arb, fmpq, fmpq_poly

root = Path(__file__).parent
widths = run_path(str(root/"ch39-current-width-normalization.py"))
inverse, study, q = widths["inverse"], widths["study"], widths["q"]


def sample(x):
    g, u, w, _, _, active = study["ports"][:6]
    return 1+(x/g**2+4/(1/(2*w**2)+active(-x)))/(18*u**2)


def enclose_shape():
    # The readings are ratios to the zero-momentum inverse response.
    nodes = [fmpq(i*i, 48**2*100) for i in range(49)]
    error = arb(0,q("1e-24").upper())
    readings = [arb(1), *(sample(arb(x))*(1+error) for x in nodes[1:])]
    polynomial = prod(fmpq_poly([-x,1]) for x in nodes)
    derivative = polynomial.derivative()
    coefficients = [arb(0) for _ in range(4)]
    for x, value in zip(nodes, readings, strict=True):
        weights = (polynomial//fmpq_poly([-x,1]))/derivative(x)
        for j in range(1,5):
            coefficients[j-1] += arb(weights[j])*value
    # For each partial fraction x/(x+r), interpolation error is exactly
    # (-1)^n r P(x)/((x+r) prod_i(r+x_i)). Use the measured end value
    # to bound the unknown sum of weights and the pole gap to bound r.
    gap, end = fmpq(1,512), nodes[-1]
    denominator = prod(gap+x for x in nodes[:-1])
    budget = (readings[-1]-1)/arb(end)
    for j in range(1,5):
        bound = sum(abs(polynomial[k])/gap**(j-k) for k in range(1,j+1))/denominator
        coefficients[j-1] += arb(0,(budget*arb(bound)).upper())
    return coefficients


def verify():
    retuned = study["retuned"]
    pole, _, (alpha,beta,family), _ = retuned["current"]["light_current"](retuned["logs"])
    assert pole > q("1/512")
    geometry = inverse["port_moments"]()
    reference = inverse["reference_coefficients"](geometry)
    estimates = enclose_shape()
    errors = []
    for value, coefficient in zip(estimates,reference[1:],strict=True):
        exact = coefficient/reference[0]
        assert value.overlaps(exact)
        errors.append(abs(value/exact-1).upper())
    print("relative coefficient errors", *(e.str(9) for e in errors))
    assert max(errors) < q("1e-9")
    width = widths["reference_width"]((alpha,beta,family))
    width *= 1+arb(0,q("1e-5").upper())
    target,L,rho,eta,c,z,_ = widths["recover"](estimates,width,geometry)
    time = study["nominal"][2]
    assert study["comp"]["model"]["target"](11430,time)>target>study["comp"]["model"]["target"](11432,time)
    masses = widths["four"]["physical_masses"](c,rho,z)
    thresholds = run_path(str(root/"ch39-neutral-current-threshold-inverse.py"))
    thresholds["exclude_lower_tops"](*masses,c,z/rho,eta,1/rho,time+L.log()/2)
    print("49 finite nodes, target interval",target.lower().str(15),target.upper().str(15))


if __name__ == "__main__":
    verify()
