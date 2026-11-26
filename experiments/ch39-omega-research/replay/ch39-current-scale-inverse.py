"""Determine the absolute scale ambiguity and invert the dimensionful current response."""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_series, fmpq_series

four = run_path(str(Path(__file__).with_name("ch39-four-current-zero-index.py")))
study, neutral, q = four["study"], four["neutral"], four["q"]
roots = study["enclose_roots"]


def verify_common_calibration(targets):
    unit, common, colour, right = neutral["graph_constants"]()
    a = [3*colour(-s, unit)/8 + right(-s, unit)/4 for s in four["momenta"]]
    readings = four["reference_readings"]()
    four["verify_precision_limit"](readings, a, common, targets)
    for eps, bounds in (("0", ("54.33050", "100.90359")),
                        ("2e-8", ("24.30561", "70.87473"))):
        ends = []
        for sign in (1, -1):
            box_corner = [v*(1+sign*(-1)**i*q(eps)) for i, v in enumerate(readings)]
            ends.append(four["recover"](box_corner, a, common)[0])
        down = (ends[0]/targets[2]).log()/2
        up = (targets[0]/ends[1]).log()/2
        for shift, floor in zip((down, up), bounds, strict=True):
            assert q(floor) < shift*10**6 < q(floor)+q(".00001")
        print("relative reading error", eps, "down/up log-scale thresholds (ppm)",
              (down*10**6).str(13), (up*10**6).str(13))


def port_moments():
    unit, common, colour, right = neutral["graph_constants"](port_builder=neutral["pole_polynomials"])
    p, r = colour
    t, v = right
    numerator, denominator = 3*p*v+2*t*r, 8*r*v
    divisor = numerator.gcd(denominator)
    assert (numerator//divisor).degree() == 10 and (denominator//divisor).degree() == 11
    def moments(pair):
        numerator, denominator = pair
        assert denominator[0] == 0 and denominator[1] != 0
        return (fmpq_series(numerator.coeffs(), prec=4)
                / fmpq_series(denominator.coeffs()[1:], prec=4))
    raw = 3*moments(colour)/8 + moments(right)/4
    A, B, C, D = [arb(raw[i])*(-1/unit)**i*(1 if i != 2 else -1) for i in range(4)]
    assert (A-q("3/70")).contains(0) and min(B, C, D) > 0
    return A, B, C/A, D/A, common/(48*arb.pi()**2)


def reference_coefficients(geometry):
    # Generate data by series division, independent of the inverse formulas.
    A, B, k, j, _ = geometry
    g, u, w = study["ports"][:3]
    z, rho = 1/(2*w*w*g*g), 1/g**2
    x = arb_series([0, 1], prec=5)
    denominator = arb_series([A, B+z, -A*k, A*j, 0], prec=5)
    return (18*u*u + rho*x + 4*rho*x/denominator).coeffs()


def recover_at_z(values, geometry, z):
    f0, f1, f2, f3 = values[:4]
    A, B, k, _, offset = geometry
    b = (B+z)/A
    scale2 = -f2*(b*b+k)/(f3*b)
    rho = -f2*A*scale2/(4*b)
    eta = f1-4*rho/A
    target = 15*scale2*(rho-offset)
    return target, scale2, rho, eta, f0/scale2


def recover(values, geometry):
    A, B, k, j, _ = geometry
    _, _, f2, f3, f4 = values
    moment_ratio = f2*f4/f3**2
    def residual(z):
        b = (B+z)/A
        return 1+(j*b-k*k)/(b*b+k)**2 - moment_ratio

    domain = arb(0).union(four["z_limit"])
    b = (B+domain)/A
    assert j*k-3*j*b*b+4*b*k*k > 0
    z = roots(residual, 0, four["z_limit"])
    target, scale2, rho, eta, c = recover_at_z(values, geometry, z)
    assert min(scale2, rho, eta, c) > 0
    four["physical_octet"](c, rho, z)
    return target, scale2, rho, eta, c, z


def verify_missing_fourth_order(values, geometry, neighbour):
    # f0, f1, f2, f3 are fixed. A nearby spectral label picks a new z.
    def gap(z):
        return recover_at_z(values, geometry, z)[0]-neighbour

    z = roots(gap, ".013", ".014") + arb(0, q("1e-12").upper())
    assert gap(z.lower()) < 0 < gap(z.upper())
    _, scale2, rho, eta, c = recover_at_z(values, geometry, z)
    assert min(scale2, rho, eta, c) > 0
    octet = four["physical_octet"](c, rho, z)
    assert q(".8417") < octet < q(".8419")


def verify_gain_ambiguity(geometry, targets, truth):
    _, _, rho, eta, c, z = truth
    for neighbour, bound in zip((targets[0], targets[2]), ("-189.699", "102.140"), strict=True):
        gain = rho/(neighbour/15 + geometry[-1])
        shift = gain.log()*10**6
        assert q(bound) < shift < q(bound)+q(".001")
        assert (15*(rho/gain-geometry[-1])-neighbour).contains(0)
        assert eta/gain > 0
        four["physical_octet"](c/gain, rho/gain, z)
        print("same full response with inverse-current log gain (ppm)", shift.str(12))


def verify():
    time = study["nominal"][2]
    targets = [study["comp"]["model"]["target"](n, time) for n in (11430,11431,11432)]
    verify_common_calibration(targets)
    geometry = port_moments()
    values = reference_coefficients(geometry)
    g, u, w = study["ports"][:3]
    truth = (targets[1], q(1), 1/g**2, 1/g**2, 18*u**2, 1/(2*w*w*g*g))
    for a, b in zip(recover(values, geometry), truth, strict=True):
        assert a.overlaps(b) and abs(a-b) < q("1e-12")
    verify_missing_fourth_order(values, geometry, targets[2])
    verify_gain_ambiguity(geometry, targets, truth)
    error = arb(0, q("1e-8").upper())
    inferred = recover([v*(1+error) for v in values], geometry)
    target, scale2, rho, eta, c, z = inferred
    assert targets[0] > target > targets[2]
    triplet, octet = four["physical_masses"](c, rho, z)
    ell = scale2.log()/2
    assert abs(ell) < q(".000002299")
    for mass, reference, bound in zip((triplet,octet), study["nominal"][:2],
                                      (".000882",".000238"), strict=True):
        assert abs(time+ell+mass.log()-reference) < q(bound)
    thresholds = run_path(str(Path(__file__).with_name("ch39-neutral-current-threshold-inverse.py")))
    thresholds["exclude_lower_tops"](triplet,octet,c,z/rho,eta,1/rho,time+ell)
    print("reference inverse-current coefficients", *(x.str(12) for x in values))
    print("inferred target limits", target.lower().str(15), target.upper().str(15))
    print("absolute log-scale error", (abs(scale2.log())/2).upper().str(10))


if __name__ == "__main__":
    verify()
