"""保留混合后的前两座族流极点，复算相邻零点的两项欧氏读数。"""
from math import comb
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, arb_poly

shape = run_path(str(Path(__file__).with_name("ch39-heavy-family-shape.py")))
mat, poles, q = (shape[key] for key in ("mat", "poles", "q"))
unit = mat["unit"]
t = q(".15")
strengths = ((".00307353742", ".00315070923"),
             (".00307528178", ".00315246483"),
             (".00308285221", ".00316008409"))


def retuned_shape(energies, z, k, h, v):
    e = energies[0]
    values = []
    for energy in energies:
        c = poles["colour"](energy, unit)/24
        r = poles["right"](energy, unit)/4
        a = z+9*c
        f = (t+e)*(t-energy)/(energy+e)+(h-k/24)*v**2/4
        g = k*energy*v**2*r
        denominator = f*(a+r)-g*a
        assert not denominator.contains(0)
        y = 4*(g-f)/denominator
        # The undressed graph has a finite, nonzero XR exchange here.
        assert not ((a+r)*y+4).contains(0) and not r.contains(0)
        values.append(y)
    d = (values[0]-values[1])/(energies[1]-e)
    b = values[0]+d*e
    assert b > 0 and d > 0
    return b, d


def check_pole_count(ports, energies, z, b, d, k, h, v):
    # Characteristic polynomial on the subspace seen by the site currents,
    # in the graph variable x = E/unit. Remove the massless root explicitly.
    (pc, vc), (pr, vr) = ([arb_poly(p) for p in pair] for pair in ports)
    energy = arb_poly([0, unit])
    e, second = energies
    y = b-energy*d
    a = vc*z-pc*(9/(24*unit))
    r = -pr*(1/(4*unit))
    f = (t-energy)*(t+e)+(energy+e)*((h-k/24)*v**2/4)
    det = f*((a*vr+r*vc)*y+4*vc*vr)-energy*(energy+e)*(k*v**2)*r*(a*y+4*vc)
    assert det[0] == 0
    poly = arb_poly(det.coeffs()[1:])
    cap = q(".11")
    assert 0 < e < second < cap < t
    for root in energies:
        assert poly(root/unit).contains(0)
    # Descartes' bound after x = (cap/unit)*u/(1+u) counts all
    # positive roots below cap, including multiplicities.
    degree = poly.degree()
    scaled = [poly[i]*(cap/unit)**i for i in range(degree+1)]
    coeffs = [sum(scaled[i]*comb(degree-i, j-i) for i in range(j+1))
              for j in range(degree+1)]
    assert all(not c.contains(0) for c in coeffs)
    signs = [c > 0 for c in coeffs]
    assert sum(a != b for a, b in zip(signs, signs[1:])) == 2


def readings(e, z, b, d, triplet, k, h, v):
    origin = mat["massive_exchange"](e, z, b, d)[2, 2]*k
    assert origin > 0
    result = []
    for log_x in (q("1/2"), q(1)):
        graph, polarization, x = shape["probe"](log_x, e, z, b, d, triplet)
        polarization[1, 1] += v**2*e*(x-1)/((t+x*e)*(t+e))
        inverse = (graph*k+mat["hypercharge"]*h).inv()-polarization
        assert mat["positive"](inverse)
        result.append(inverse.inv()[2, 2]/(x*origin))
    return result


def verify():
    time = poles["study"]["nominal"][2]
    target = poles["study"]["comp"]["model"]["target"](11430, time)
    *_, colour, right = poles["neutral"]["graph_constants"](
        port_builder=poles["neutral"]["pole_polynomials"])
    for case, logs in enumerate(poles["retune"]["log_sets"]):
        model, _, _ = poles["retune"]["predictions"](logs, time)
        scale, z, rho, eta, c, masses = model
        b0, d0 = c/rho, eta/rho
        n1 = poles["reference_readings"](logs)[3]
        n2 = mat["geometry"]["second_neutral_mass"](scale, z, b0, d0)
        e, second = n1**2/scale, n2**2/scale
        h_ref = mat["flow"]["running"]["hypercharge_window"](time+n1.log()).mid()**2
        refs = readings(e, z, b0, d0, masses[0], 1/rho, h_ref, q(0))
        k = 1/(poles["offset"]+target/(15*scale))
        limits = (("48.42", "16.36"), ("48.49", "16.49"), ("48.60", "16.66"))
        for strength, h, limit in zip(strengths[case], (h_ref, q(".01")), limits[case], strict=True):
            v = q(strength)
            b, d = retuned_shape((e, second), z, k, h, v)
            check_pole_count((colour, right), (e, second), z, b, d, k, h, v)
            kinetic = arb_mat([[k/24, 0, 0], [0, k/4, 0], [0, 0, k/d]])
            kinetic += mat["hypercharge"]*(h-k/24)
            cost = v**2*kinetic[1, 1]/(t+e)
            assert mat["positive"](kinetic) and 0 < cost < 1
            triplet = masses[0]*(b/b0).sqrt()
            assert 4*triplet**2 < t
            mass_shift = (triplet/masses[0]-1)*10**6
            assert q("-18") < mass_shift < q("-16.9")
            amps = readings(e, z, b, d, triplet, k, h, v)
            errors = [(value/ref-1)*10**9 for value, ref in zip(amps, refs, strict=True)]
            assert -q(limit) < errors[0] < 0 < errors[1] < q(limit)
            print("case", case, "g_Y", h.sqrt().str(10), "v_R", strength,
                  "B/B0", (b/b0).str(12), "D/D0", (d/d0).str(12),
                  "kinetic fraction", cost.str(8), "triplet shift (ppm)", mass_shift.str(8), "errors (ppb)",
                  *(error.str(10) for error in errors), flush=True)


if __name__ == "__main__":
    verify()
