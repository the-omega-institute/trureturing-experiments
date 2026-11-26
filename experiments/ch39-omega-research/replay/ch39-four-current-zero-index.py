"""Identify the arithmetic target from four samples of the neutral exchange kernel."""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import arb

study = run_path(str(Path(__file__).with_name("ch39-top-threshold-resolution.py")))
neutral = study["retuned"]["current"]["neutral"]
q = study["q"]
momenta = tuple(map(q, ("1/1024", "1/320", "1/160")))
z_limit = (-q("2/21")).exp() / 56


def reference_readings():
    g, u, w, _, _, active = study["ports"][:6]
    c = 18 * u**2
    beta = 1 / (2 * w**2)
    return [1/c, *(1/(c + s/g**2 + 4/(beta + active(-s))) for s in momenta)]


def recover(readings, coefficients, common):
    s1, s2, s3 = momenta
    a1, a2, a3 = coefficients
    alpha12 = (s2*a2 - s1*a1) / (s2 - s1)
    alpha23 = (s3*a3 - s2*a2) / (s3 - s2)
    assert min(a1, a2, a3, alpha12, alpha23) > 0
    factor = (s2 - s1)*s3 / ((s3 - s2)*s1)

    def ratio(z):
        return factor*(z + alpha12)*(z + a3) / ((z + alpha23)*(z + a1))

    domain = arb(0).union(z_limit)
    slope = ((alpha23-alpha12)/((domain+alpha12)*(domain+alpha23))
             + (a1-a3)/((domain+a3)*(domain+a1)))
    assert slope > 0
    c, *positive = [1/k for k in readings]
    v1, v2, v3 = [(f-c)/s for f, s in zip(positive, momenta, strict=True)]
    assert v1 > v2 > v3
    observed = (v1-v2)/(v2-v3)
    z = study["enclose_roots"](lambda value: ratio(value)-observed, 0, z_limit)
    # Factor the difference to avoid subtracting two overlapping intervals.
    d1, d3 = s1*(z+a1), s3*(z+a3)
    delta = (s3-s1)*z + s3*a3-s1*a1
    rho = (v1-v3)*d1*d3/(4*delta)
    eta = v3 - 4*rho/d3
    assert min(rho, eta) > 0
    target = 15*(rho - common/(48*arb.pi()**2))
    return target, rho, eta, z


def verify_three_point_witness(readings, coefficients, rho, indices):
    """At fixed alternative site coupling, fit the three retained inverse readings."""
    s = (arb(0), *momenta)
    a = (None, *coefficients)
    observed = [1/k for k in readings]
    i, j, k = indices

    def line(z):
        def y(index):
            return observed[index] - (4*rho/(z+a[index]) if index else 0)
        eta = (y(k)-y(i))/(s[k]-s[i])
        c = y(i)-eta*s[i]
        return c+eta*s[j]-y(j), c, eta

    z = study["enclose_roots"](lambda x: line(x)[0], ".013", z_limit)
    # Prove existence by endpoint signs even without invoking root uniqueness.
    z += arb(0, q("1e-10").upper())
    assert 0 < z < z_limit
    assert line(z.lower())[0] < 0 < line(z.upper())[0]
    _, c, eta = line(z)
    assert min(c, eta) > 0
    return physical_octet(c, rho, z)


def physical_masses(c, rho, z):
    d = ((1 + q("21/2")*(56*z).log())/5).exp()
    def mass_equation(h):
        a_squared = ((1+4*d.log()+h.log()+arb(15).log()/2)/10).exp()
        return 15*a_squared - 1/(4*z) - (d/h)**2

    h = study["enclose_roots"](mass_equation, ".05", "64")
    triplet = h*(c/(18*rho)).sqrt()
    octet = (triplet**2+d**2).sqrt()
    assert 0 < triplet < octet < 1
    return triplet, octet


def physical_octet(c, rho, z):
    return physical_masses(c, rho, z)[1]


def verify_precision_limit(readings, coefficients, common, neighbours):
    # On the whole input box there is one smooth physical inverse.
    error = arb(0, q("7e-8").upper())
    box = [k*(1+error) for k in readings]
    _, rho, _, z = recover(box, coefficients, common)
    physical_octet(1/box[0], rho, z)
    # Cofactors for rho remove its column from the Jacobian of F. Each
    # remaining determinant is rho times a chord test for T'(z).
    s = (arb(0), *momenta)
    slope = [arb(0), *[-4/(z+a)**2 for a in coefficients]]
    for i, j, k in combinations(range(4), 3):
        minor = (s[j]-s[i])*(slope[k]-slope[i]) - (s[k]-s[i])*(slope[j]-slope[i])
        assert minor < 0
    # The four signs alternate. The two opposite alternating corners
    # therefore attain both extremes, with orientation certified below.
    for neighbour, direction, lower, upper in (
        (neighbours[2], 1, "36.19224668", "36.19224670"),
        (neighbours[0], -1, "67.19390529", "67.19390532"),
    ):
        def crossing(ppb):
            values = [k*(1+direction*(-1)**i*ppb*q("1e-9")) for i, k in enumerate(readings)]
            target = recover(values, coefficients, common)[0]
            return direction*(neighbour-target)

        assert crossing(q(lower)) < 0 < crossing(q(upper))
        print("critical relative error in ppb between", lower, "and", upper)


def verify():
    unit, common, colour, right = neutral["graph_constants"]()
    coefficients = [3*colour(-s, unit)/8 + right(-s, unit)/4 for s in momenta]
    readings = reference_readings()
    target_at = study["comp"]["model"]["target"]
    time = study["nominal"][2]
    neighbours = [target_at(n, time) for n in (11430, 11431, 11432)]
    nominal, g, w = neighbours[1], study["ports"][0], study["ports"][2]
    exact = recover(readings, coefficients, common)
    truth = (nominal, 1/g**2, 1/g**2, 1/(2*w**2*g**2))
    for value, expected in zip(exact, truth, strict=True):
        assert value.overlaps(expected) and abs(value-expected) < q("1e-12")
    error = arb(0, q("2e-8").upper())
    inferred, rho, eta, z = recover([k*(1+error) for k in readings], coefficients, common)
    assert neighbours[0] > q("3.29442601") > inferred > q("3.29372669") > neighbours[2]
    for neighbour in (neighbours[0], neighbours[2]):
        alternative_rho = neighbour/15 + common/(48*arb.pi()**2)
        for indices in combinations(range(4), 3):
            octet = verify_three_point_witness(readings, coefficients, alternative_rho, indices)
            assert q(".8764") < octet < q(".9141")
    verify_precision_limit(readings, coefficients, common, neighbours)
    print("four readings", *(k.str(17) for k in readings))
    print("neighbour targets", *(t.str(15) for t in neighbours))
    print("recovered target", inferred.lower().str(16), inferred.upper().str(16))
    print("log errors in g, g_X, z", *(abs((v/v0).log()/p).upper().str(8)
          for v, v0, p in zip((rho, eta, z), truth[1:], (2, 2, 1), strict=True)))


if __name__ == "__main__":
    verify()
