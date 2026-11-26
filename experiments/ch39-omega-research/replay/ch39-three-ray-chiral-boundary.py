"""同一组三条有理射线上的手征边界、内部零模及重层有限阈值。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx

arc = run_path(str(Path(__file__).with_name("ch39-primitive-pole-arc-readout.py")))
model = arc["model"]
soft = model["shell"]["soft"]
geo = soft["geo"]


def boundary_sites():
    sites = geo["sites"]
    boundary = []
    for axis in range(3):
        a, b = [i for i in range(3) if i != axis]
        # On the positive ray from the pole towards (5,1,1)/3,
        # x_axis increases relative to the other two equal coordinates.
        ray = [i for i, v in enumerate(sites) if v[a] == v[b] and v[a] < v[axis] <= 5*v[a]]
        radii = {i: sum((x-17)**2 for x in sites[i]) for i in ray}
        closest = [i for i in ray if radii[i] == min(radii.values())]
        assert len(closest) == 1 and radii[closest[0]] == 68
        boundary.extend(closest)
    assert len(set(boundary)) == 3
    return boundary


def laplacian(charge, height, *, phases=None):
    matrix = acb_mat(19, 19)
    if phases is None:
        phases = geo["phases"](height)
    pole, tuned = soft["pole"], soft["star"]["orbit"]
    delta = arb(23)/25*arb(2)**-20/(height**2 + soft["point"]["study"]["offset"])
    for i, j, d, z in geo["edges"]:
        spoke = pole in (i, j)
        weight = 1-delta if spoke and (i in tuned or j in tuned) else 1
        phase = acb(1) if spoke else (z*phases.get(d, 1))**charge
        matrix[i,i] += weight
        matrix[j,j] += weight
        matrix[i,j] -= weight*phase
        matrix[j,i] -= weight*phase.conjugate()
    return matrix


def high_squares(matrix, boundary):
    mass_block = acb_mat([[matrix[i,j] for j in range(19)] for i in range(19) if i not in boundary])
    squared = mass_block*mass_block.conjugate().transpose()
    levels = [v.real for v in squared.eig(multiple=True)]
    assert len(levels) == 19-len(boundary) and all(x > 0 for x in levels)
    return levels


def squared_mass_bounds(flat, boundary):
    keep = [i for i in range(19) if i not in boundary]
    n = len(keep)
    # Every interior vertex is joined to every boundary vertex. The
    # principal magnetic Laplacian H thus dominates this diagonal floor.
    floor = min(sum(abs(flat[i,k]) for k in boundary).lower() for i in keep)
    assert floor > arb(len(boundary))*(1-arb("1e-12"))
    mean = sum(flat[i,i].real for i in keep)/n
    scatter = sum((flat[i,i].real-mean)**2 for i in keep)
    scatter += sum(abs(flat[i,j])**2 for i in keep for j in keep if i != j)
    border = sum(abs(flat[i,k])**2 for i in keep for k in boundary)
    # Trace and squared deviation are phase independent. Cauchy-Schwarz on
    # the two groups around each ordered eigenvalue bounds it on both sides.
    # D D* = H² + B B*, and ||B B*|| <= trace(B B*).
    bounds = []
    for j in range(n):
        low = (mean-(scatter*(n-1-j)/(n*(j+1))).sqrt()).lower()
        high = min(arb(36), (mean+(scatter*j/(n*(n-j))).sqrt()).upper())
        bounds.append((max(floor,low)**2).union(high**2 + border))
    return bounds


def uniform_exclusion(height, boundary, pole, vectors, radius, u):
    flat = laplacian(0, height)
    fermions = squared_mass_bounds(flat, boundary)
    scalars = squared_mass_bounds(flat, pole)
    # Signs of P times each heavy coefficient come directly from the group
    # indices. Increasing a squared mass decreases its finite exposure.
    projected = lambda row: -(115*row[0]-333*row[1]+218*row[2])/218
    for charge,row in zip(model["charges"],model["weyl"]):
        sign = 1 if charge in (1,3) else -1
        assert sign*projected(row) > 0
    assert projected(model["higgs"]) > 0
    inverses = []
    for maximize_inverse in (False, True):
        levels = {0:vectors}
        for q in model["charges"]:
            heavy = (q in (1,3)) == maximize_inverse
            levels[q] = [x.upper() if heavy else x.lower() for x in fermions]
        higgs = [x.upper() if maximize_inverse else x.lower() for x in scalars]
        inverse, _, _ = arc["readout"](levels, radius, u, fermion_weight=2,
                                       higgs_levels=higgs, min_mu="1e10")
        inverses.append(inverse)
    strong = 1/inverses[0].union(inverses[1])
    assert arb(".07058") < strong < arb(".08201")
    print("任意固定物质相位的强耦合必要范围", 1/inverses[1], 1/inverses[0])


def verify():
    ctx.prec = 224
    height = acb.zeta_zero(11431).imag
    boundary = boundary_sites()
    pole = [soft["pole"]]
    levels = {q: high_squares(laplacian(q,height), boundary) for q in model["charges"]}
    levels[0] = model["shell"]["shell_spectra"]((0,))[0]
    higgs = high_squares(laplacian(3,height), pole)
    radius = geo["theta"](804)/geo["theta"](12)

    def read(u, min_mu):
        return arc["readout"](levels, radius, u, fermion_weight=2, higgs_levels=higgs, min_mu=min_mu)

    geometric = arc["first_arc"]()
    uniform_exclusion(height, boundary, pole, levels[0], radius, geometric)
    inv, mu, common = read(geometric, "1e10")
    assert arb(".07442382") < 1/inv < arb(".07442383")
    assert arb("1.58746e12") < mu < arb("1.58747e12") and common > 0
    print("几何刚度下的强耦合、电弱交点 GeV、共同电弱逆耦合", 1/inv, mu, common)
    def projection(u):
        # Root for the book's fixed central alpha_s.
        return 1/arb(".1180") - read(u, "1e13")[0]
    lower, upper = arb(1), arb(4096)
    assert projection(lower)<0<projection(upper)
    for _ in range(40):
        middle = (lower+upper)/2
        if projection(middle)<0:
            lower = middle
        else:
            assert projection(middle)>0
            upper = middle
    u = lower.union(upper)
    _, scale, coupling = read(u, "1e13")
    assert arb("51.7005") < u < arb("51.7006")
    assert arb("5.1028e14") < scale < arb("5.1030e14") and coupling > 0
    print("中心三支匹配所需刚度、尺度 GeV、共同逆耦合", u, scale, coupling)


if __name__ == "__main__":
    verify()
