"""用两个整数荷矢量型 Dirac 场检验有限动量族流的零点可辨识性。"""
from math import comb
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat

shape = run_path(str(Path(__file__).with_name("ch39-heavy-family-shape.py")))
mat, poles, q = (shape[key] for key in ("mat", "poles", "q"))
mixed_charge = arb_mat([[6], [7], [-4]])
centres = ((".00562963", ".25208230"), (".00559847", ".25217338"),
           (".00546305", ".25257022"))
three_point = ((".00430813", ".26449002"), (".00427503", ".26459128"),
               (".00413118", ".26503244"))


def dirac(x, r):
    """一次减法的单位荷极化和对 r 的导数，含交错级数余项。"""
    assert x > 1 and r >= 0 and x*r < 1
    value = slope = arb(0)
    for n in range(1, 29):
        coefficient = (x**n-1)/(n*(2*n+3)*comb(2*n+2, n+1))
        term = coefficient*r**(n-1)
        signed = term if n % 2 else -term
        value += signed*r
        slope += n*signed
    n = 29
    next_term = (x**n-1)*r**(n-1)/(n*(2*n+3)*comb(2*n+2, n+1))
    factor = 2*arb.pi()**2
    return ((value+arb(0).union(next_term*r))/factor,
            (slope+arb(0).union(n*next_term))/factor)


def packet(log_x, e, z, b, d, triplet, k0, k, h):
    p, a, x = shape["probe"](log_x, e, z, b, d, triplet)
    hyper = mat["hypercharge"]*h
    bare = p*k+hyper
    inverse = bare.inv()-a
    assert mat["positive"](inverse)
    prop = inverse.inv()
    jcol = arb_mat([[p[i, 2]] for i in range(3)])
    w0, w = ((a.inv()-p*ki-hyper).inv() for ki in (k0, k))
    # For the normalized family entry, the tree part cancels exactly.
    gap = (jcol.transpose()*w0*(jcol+p*w*jcol*k))[0, 0]*k*(k0-k)
    base = prop[2, 2]
    target = base+gap
    independent = ((p*k0+hyper).inv()-a).inv()[2, 2]*(k/k0)
    assert target.overlaps(independent) and gap > 0
    v = prop*mixed_charge
    return x, base, target, gap, v[2, 0]**2, (mixed_charge.transpose()*v)[0, 0]


def equation(data, r1, r2):
    x, base, target, gap, overlap, self_entry = data
    f1, slope1 = dirac(x, r1)
    f2, slope2 = dirac(x, r2)
    denominator = 1-f2*self_entry
    assert denominator > 0
    gain = f2*overlap/denominator
    enlarged = base+gain
    last = 1-f1*enlarged
    assert last > 0
    residual = (f1*target*enlarged+gain-gap)/target
    gradient = [slope1*enlarged, slope2*overlap*(1/target+f1)/denominator**2]
    return residual, gradient, residual/last


def certify(packets, centre, inverse=None):
    # A fixed rational matrix is sufficient for all three contraction tests.
    if inverse is None:
        inverse = arb_mat([[q("-.074"), q(".032")], [q(".760"), q("-.232")]])*10**6
    radius = q(".000001")
    error = arb(0, radius.upper())
    box = [value+error for value in centre]
    residual = arb_mat([[equation(row, *centre)[0]] for row in packets])
    gradient = arb_mat([equation(row, *box)[1] for row in packets])
    defect = arb_mat([[1, 0], [0, 1]])-inverse*gradient
    norm = max(sum(abs(defect[i, j]) for j in range(2)).upper() for i in range(2))
    offset = max(abs(value).upper() for value in (inverse*residual).entries())
    assert norm < q(".061") and offset+norm*radius < radius and not inverse.det().contains(0)
    return box, norm


def hypercharge_response(radii, params):
    e, z, b, d, triplet, k0, k, h = params
    p, a, x = shape["probe"](q(1), e, z, b, d, triplet)
    start = mat["massive_exchange"](e, z, b, d)
    y = mat["y"]
    f1, f2 = [dirac(x, r)[0] for r in radii]
    extra = mixed_charge*mixed_charge.transpose()*f2
    extra[2, 2] += f1
    ratios = []
    for coupling, vacuum in ((k0, a), (k, a+extra)):
        inverse = (p*coupling+mat["hypercharge"]*h).inv()-vacuum
        assert mat["positive"](inverse)
        origin = start*coupling+mat["hypercharge"]*h
        ratios.append((y.transpose()*inverse.inv()*y)[0, 0]/(x*(y.transpose()*origin*y)[0, 0]))
    return ratios[1]/ratios[0]-1


def all_mass_bound(measurements, params):
    # For a given second mass, the first fermion contributes only to the XX
    # element of the inverse. Invert the target interval in this scalar channel.
    limit = q(".33")
    e, _, _, _, triplet, *_ = params
    assert e/triplet**2 < limit
    cells = 256
    grid = [limit*i/cells for i in range(cells+1)]
    tolerance = q("2e-8")
    tables = [[dirac(data[0], r)[0] for r in grid] for data in measurements]

    def needed(data, f2, error):
        _, base, target, gap, overlap, self_entry = data
        den = 1-f2*self_entry
        assert den > 0
        increment = f2*overlap/den
        return (gap+target*error-increment)/((1+error)*target*(base+increment))

    candidates = []
    for j in range(cells):
        windows = [(needed(data, table[j+1], -tolerance), needed(data, table[j], tolerance))
                   for data, table in zip(measurements, tables, strict=True)]
        for i in range(cells):
            if all(not (table[i+1] < low or table[i] > high)
                   for table, (low, high) in zip(tables, windows, strict=True)):
                candidates.append(j)
                break
    assert candidates
    lower_r2 = grid[min(candidates)]
    # Positivity at the upper corner establishes the monotone positive
    # propagator branch on the entire square. Test the two family points too.
    hypercharge_response((limit, limit), params)
    for data in measurements:
        equation(data, limit, limit)
    lower_y = hypercharge_response((q(0), lower_r2), params)
    print("all masses: lower r2", lower_r2.str(9), "hypercharge lower bound", lower_y.lower().str(9))
    return lower_y


def verify():
    time = poles["study"]["nominal"][2]
    goal = poles["study"]["comp"]["model"]["target"](11430, time)
    assert (mat["y"].transpose()*mixed_charge)[0, 0].contains(q("4.5"))
    for case, logs in enumerate(poles["retune"]["log_sets"]):
        (scale, z, rho, eta, c, masses), _, _ = poles["retune"]["predictions"](logs, time)
        charged, _, _, neutral, _ = poles["reference_readings"](logs)
        e = neutral**2/scale
        triplet = masses[0]
        h = mat["flow"]["running"]["hypercharge_window"](time+neutral.log()).mid()**2
        k = 1/(poles["offset"]+goal/(15*scale))
        params = e, z, c/rho, eta/rho, triplet, 1/rho, k, h
        measurements = [packet(t, *params) for t in (q("1/2"), q(1))]
        box, norm = certify(measurements, [q(v) for v in centres[case]])
        mass_ratios = [(e/r).sqrt()/triplet for r in box]
        assert q("7.58") < mass_ratios[0] < q("7.71") and q("1.13") < mass_ratios[1] < q("1.14")
        low = packet(q("1/4"), *params)
        third = equation(low, *box)[2]*10**9
        assert q("-20.2") < third < q("-20.03")
        radii = [q(v) for v in three_point[case]]
        assert all((e/v).sqrt() > triplet for v in radii)
        shifts = [equation(row, *radii)[2]*10**9 for row in (low, *measurements)]
        assert all(abs(s) < q("12.53") for s in shifts)
        y_shift = q("4.5")**2*dirac(measurements[-1][0], box[1])[0]
        assert q(".0134") < y_shift < q(".0136")
        running_shift = -27*(charged/(scale*e/box[1]).sqrt()).log()/(8*arb.pi()**2)
        assert (scale*e/box[1]).sqrt() < charged < (scale*e/box[0]).sqrt()
        assert q("-.206") < running_shift < q("-.204")
        comparisons = [hypercharge_response(uv, params) for uv in (box, radii)]
        assert all(q(".00204") < response < q(".00214") for response in comparisons)
        global_lower = all_mass_bound(measurements, params)
        assert global_lower > q(".00167")
        print("case", case, "mass ratios", *(value.str(10) for value in mass_ratios),
              "contraction", norm.str(5), "third relative shift (ppb)", third.str(9),
              "three readouts (ppb)", *(s.str(8) for s in shifts),
              "hypercharge responses", *(s.str(8) for s in comparisons),
              "additional hypercharge polarization", y_shift.str(9),
              "step running", running_shift.str(9), flush=True)


if __name__ == "__main__":
    verify()
