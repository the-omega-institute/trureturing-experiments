"""两个夸克跃迁概率的加权和在参考六面相位处严格极小。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, arb_poly, arb_series, ctx


root = Path(__file__).parent
invariant = run_path(str(root / "ch39-flavor-invariants.py"))
joint = run_path(str(root / "ch39-common-vub-weights.py"))


def spectral_data():
    flavor = invariant["flavor"]
    matrix = flavor["reference_matrix"]()
    unit = arb(".01401")**2
    characteristic = [p*(1/unit**i) for i, p in enumerate(invariant["invariants"]()[:3], 1)]
    nodes = [arb(j) for j in (-1, 0, 1)]
    data = []
    for c in nodes:
        y = matrix(invariant["phase"](c))
        h = (y*y.conjugate().transpose())*(1/unit)
        squared = h*h
        data.append([h[0,0].real, h[1,1].real, squared[0,0].real, squared[1,1].real])
    coefficients = arb_mat([[1, v, v*v] for v in nodes]).solve(arb_mat(data))
    diagonals = [arb_poly([coefficients[i,j] for i in range(2 if j < 2 else 3)]) for j in range(4)]
    return characteristic + diagonals


def evaluate(poly, x):
    value = 0
    for coefficient in reversed(poly.coeffs()):
        value = value*x + coefficient
    return value


def probabilities(c, data):
    x = arb_series([c, 1], prec=3)
    e1, e2, e3, h0, h1, s0, s1 = [evaluate(poly, x) for poly in data]

    def cubic(t):
        return ((t-e1[0])*t+e2[0])*t-e3[0]

    def derivative(t):
        return (3*t-2*e1[0])*t+e2[0]

    low, high = arb(".85"), arb("1.15")
    bounds = low.union(high)
    assert cubic(low) < 0 < cubic(high) and derivative(bounds) > 0 and 6*low-2*e1[0] > 0
    for _ in range(5):
        mid = bounds.mid()
        bounds = bounds.intersection(mid-cubic(mid)/derivative(bounds))
    linear = (e1[1]*bounds**2-e2[1]*bounds+e3[1])/derivative(bounds)
    t = arb_series([bounds, linear], prec=3)
    quadratic = -(((t-e1)*t+e2)*t-e3)[2]/derivative(bounds)
    t = arb_series([bounds, linear, quadratic], prec=3)
    denominator = (3*t-2*e1)*t+e2
    return [(h*(t-e1)+s+e3/t)/denominator for h, s in ((h0, s0), (h1, s1))]


def verify():
    ctx.prec, ctx.cap = 320, 3
    model = joint["model"]
    flavor = invariant["flavor"]
    reference = flavor["singlet"](11661)[1]
    data = spectral_data()
    up, charm = probabilities(reference, data)
    ratio = -up[1]/charm[1]
    assert arb("1.16626614") < ratio < arb("1.16626615")
    count = 512
    for i in range(count):
        left, right = [arb(-1)+2*arb(j)/count for j in (i, i+1)]
        u, v = probabilities(left.union(right), data)
        half = u[2]+ratio*v[2]
        assert arb("3.5297e-8") < half < arb("3.5365e-8")
    print("两个跃迁概率的支撑直线系数", ratio.str(15))

    shared = model["shared"]
    times = shared["clock"]["clock"]["times"]()
    _, _, lambdas, slopes = shared["flavor_data"]()
    columns, charm_readings, squares = [], [], []
    read = flavor["reference_readout"]()
    for n in joint["indices"]:
        height = acb.zeta_zero(n).imag
        phase, c = flavor["singlet"](n)
        _, s23, s13, _ = read(phase)
        pu, pc = probabilities(c, data)
        u, v = s13**2, s23**2*(1-s13**2)
        assert pu[0].overlaps(u) and pc[0].overlaps(v)
        x = (joint["scale"]/height)**2
        columns.append([arb(1), *[1/(1+j*j*x) for j in range(1, 5)], model["gradient"](x, lambdas, slopes), u])
        squares.append((c-reference)**2)
        charm_readings.append(v)
    weights = arb_mat(list(zip(*columns, strict=True))).solve(arb_mat([[v] for v in [*times[:5], 0, up[0]]])).entries()
    assert all(w > 0 for w in weights)
    variance = sum(w*s for w, s in zip(weights, squares, strict=True))
    excess = sum(w*v for w, v in zip(weights, charm_readings, strict=True))-charm[0]
    assert arb("5.67e-13") < excess < arb("5.70e-13")
    assert excess.overlaps(arb("3.5297e-8").union(arb("3.5365e-8"))*variance/ratio)
    print("七点解的第二个夸克跃迁概率超额", excess.str(16))


if __name__ == "__main__":
    verify()
