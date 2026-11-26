"""全体非辐边共用一个刚度时，指定矢量阈值的唯一解及味基态反例。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx

star = run_path(str(Path(__file__).with_name("ch39-neutral-clock-star-threshold.py")))
geo, point, pole = star["geo"], star["point"], star["pole"]
n = len(geo["sites"])
identity = acb_mat([[int(i == j) for j in range(n)] for i in range(n)])


def vector_pencil(height):
    base, neutral, charged = (acb_mat(n, n) for _ in range(3))
    phases = geo["phases"](height)
    for i, j, d, z in geo["edges"]:
        if pole in (i, j):
            pairs = [(base, acb(1))]
        else:
            pairs = [(neutral, acb(1)), (charged, (z*phases.get(d, 1))**5)]
        for matrix, link in pairs:
            matrix[i, i] += 1
            matrix[j, j] += 1
            matrix[i, j] = -link
            matrix[j, i] = -link.conjugate()
    return base, neutral, charged


def simple_root(height):
    base, neutral, charged = vector_pencil(height)
    lifted = base+acb_mat([[arb(1)/n]*n for _ in range(n)])
    radius = geo["theta"](804)/geo["theta"](12)
    target = arb("3.664")
    factor = radius**2*(-8*arb.pi()*target/7).exp()
    def polynomial(z):
        return (base+z*charged).det()-factor*(lifted+z*neutral).det()
    # The pole row in each variable matrix vanishes, so degree is at most 18.
    samples = [polynomial(acb(0, 2*arb.pi()*k/20).exp()) for k in range(20)]
    coefficients = [v.real/20 for v in acb.dft(samples)]
    assert coefficients[0] < 0 and all(c > 40 for c in coefficients[1:19])
    lo, hi = arb("1e-8"), arb(".1")
    assert polynomial(lo).real < 0 and polynomial(hi).real > 0
    for _ in range(96):
        mid = (lo+hi)/2
        sign = polynomial(mid).real
        if sign < 0:
            lo = mid
        else:
            assert sign > 0
            hi = mid
    weight = lo.union(hi)
    return weight, base+weight*charged


def flavor_block(h, weight):
    m = geo["matrix"](0, geo["phases"](h))
    for i, j, _, _ in geo["edges"]:
        if pole not in (i, j):
            m[i, i] += weight-1
            m[j, j] += weight-1
            m[i, j] *= weight
            m[j, i] *= weight
    native = point["native"]
    diag = acb_mat([[m[k, 0] if i == j else 0 for j in range(2)]
                    for i, k in enumerate(native["portals"])])
    source = native["source"]*diag*native["columns"].transpose()*star["flavor_port"]()
    for i, v in enumerate(native["inside"]):
        m[v, 0] = star["gain"]*source[i, 0]
        m[0, v] = m[v, 0].conjugate()
    m[0, 0] = -star["gain"]**2/(h*h+point["study"]["offset"])
    return m


def lowest(matrix, lo, hi):
    def positive_at(value):
        shifted = matrix-value*identity
        for k in range(1, n+1):
            minor = acb_mat([[shifted[i, j] for j in range(k)] for i in range(k)])
            sign = minor.det().real
            if sign < 0:
                return False
            assert sign > 0
        return True
    assert positive_at(lo) and not positive_at(hi)
    for _ in range(62):
        mid = (lo+hi)/2
        if positive_at(mid):
            lo = mid
        else:
            hi = mid
    return lo.union(hi)


def verify():
    ctx.prec = 224
    h = acb.zeta_zero(11431).imag
    weight, charged = simple_root(h)
    light = lowest(charged, arb(0), arb(".01"))
    center = lowest(flavor_block(h, weight), arb(-10), arb(-8))
    competitor = lowest(flavor_block(acb.zeta_zero(1069).imag, weight), arb(-10), arb(-8))
    assert arb(".00006552972566593") < weight < arb(".00006552972566594")
    assert arb(".00014101149364790") < light < arb(".00014101149364791")
    assert center-competitor > arb("2.4683e-11")
    print("Fitted stiffness, squared mass, zero 1069 flavor advantage",
          weight.str(17), light.str(17), (center-competitor).str(12))

    # This linear combination kills both gauge index vectors (5,3,2) and (0,2,3).
    direction = (1, -3, 2)
    for indices in ((5,3,2), (0,2,3), (1,1,1)):
        assert sum(x*y for x, y in zip(direction, indices)) == 0
    em, weak, color = arb("127.930"), arb(".23129"), arb(".1180")
    couplings = [3*(1-weak)*em/5, weak*em, 1/color]
    beta = [arb(41)/10, -arb(19)/6, -7]
    initial = sum(x*y for x, y in zip(direction, couplings))
    rate = -sum(x*y for x, y in zip(direction, beta))/(2*arb.pi())
    log_range = (arb("1e13")/arb("91.1876")).log().union((arb("1e17")/arb("91.1876")).log())
    residual = initial+rate*log_range
    assert residual < -arb("10.6083") and -residual/3 > arb("3.5361")
    print("Untouched one-loop contrast over the matching window", residual.lower().str(8),
          residual.upper().str(8))


if __name__ == "__main__":
    verify()
