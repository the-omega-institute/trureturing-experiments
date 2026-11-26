"""用一个暗模与一个较重的右弱矢量构造重配中性图的精确匹配。"""
from math import comb
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, arb_poly

root = Path(__file__).parent
ordered = run_path(str(root/"ch39-ordered-neutral-mass-bound.py"))
retuned = run_path(str(root/"ch39-retuned-local-vector.py"))
poles, mat, q = (retuned[key] for key in ("poles", "mat", "q"))
unit = mat["unit"]
eye = arb_mat([[1, 0], [0, 1]])
family = arb_mat([[0], [0], [1]])
y = arb_mat([[q("1/6")], [q("1/2")]])
centres = ((".00003290611134", ".00006086317196"),
           (".00003293989732", ".00006093452819"),
           (".00003308608722", ".00006124346493"))


def ports_at(energy):
    return poles["colour"](energy, unit)/24, poles["right"](energy, unit)/4


def tensor(v):
    return v*v.transpose()


def scalar(left, matrix, right):
    return (left.transpose()*matrix*right)[0, 0]


def dark_atom(ports, h, k, triplet):
    (pc, vc), (pr, vr) = ([arb_poly(a) for a in pair] for pair in ports)
    det = arb_poly([0, 1])*pc*pr+(h/k-q("1/24"))*(vc*pr*q("2/3")+pc*vr)
    assert det[0] == 0
    poly = arb_poly(det.coeffs()[1:])
    left, right = [q(m)**2*triplet**2/unit for m in ("2.08", "2.10")]
    assert ordered["variations_in"](poly, left, right) == 1
    sign = 1 if poly(right) > 0 else -1
    x = poles["roots"](lambda v: sign*poly(v), left, right)
    mid = x.mid()
    x = x.intersection(mid-poly(mid)/poly.derivative()(x))
    t = unit*x
    c, r = ports_at(t)
    v = arb_mat([[1/(6*c)], [1/(2*r)]])
    return t, v*(1/(v.transpose()*v)[0, 0].sqrt())


def fit(energies, z, k, h, atoms, strengths):
    values = []
    slopes = []
    e, other = energies
    for E in energies:
        c, r = ports_at(E)
        direction = arb_mat([[3*c], [-r]])
        bare = arb_mat([[-k*E*c, 0], [0, -k*E*r]])+(h-k/24)*tensor(y)
        factors = [(E+e)/((t+e)*(t-E)) for t, _ in atoms]
        added = sum((a*u*tensor(v) for a, u, (_, v) in zip(factors, strengths, atoms, strict=True)),
                    arb_mat(2, 2))
        leg = (eye+bare*added).inv()*direction
        effective = z+9*c+r+k*E*scalar(direction, added, leg)
        value = -4/effective
        values.append(value)
        slopes.append([value**2*k*E*a*(v.transpose()*leg)[0, 0]**2/4
                       for a, (_, v) in zip(factors, atoms, strict=True)])
    d = (values[0]-values[1])/(other-e)
    ds = [(a-b)/(other-e) for a, b in zip(*slopes, strict=True)]
    return values[0]+e*d, d, [a+e*b for a, b in zip(slopes[0], ds, strict=True)], ds


def integral_slope(r):
    assert 0 < r < 1
    terms = [r**(n-1)/((2*n+3)*(2*n+1)*comb(2*n, n)) for n in range(1, 26)]
    return sum((-1)**j*v for j, v in enumerate(terms[:-1]))+arb(0).union(terms[-1])


def readouts(data, strengths):
    energies, z, k, h, triplet0, b0, atoms, refs = data
    e = energies[0]
    b, d, db, dd = fit(energies, z, k, h, atoms, strengths)
    triplet = triplet0*(b/b0).sqrt()
    vectors = [arb_mat([[v[0, 0]], [v[1, 0]], [0]]) for _, v in atoms]
    values, derivs = [], []
    for log_x in (q(0), q(".5"), q(1)):
        x = log_x.exp()
        s = x*e
        if log_x == 0:
            tree, A = mat["massive_exchange"](e, z, b, d), arb_mat(3, 3)
        else:
            tree, A, _ = retuned["shape"]["probe"](log_x, e, z, b, d, triplet)
        factors = [(s-e)/((t+e)*(t+s)) for t, _ in atoms]
        for u, factor, v in zip(strengths, factors, vectors, strict=True):
            A += u*factor*tensor(v)
        bare = tree*k+mat["hypercharge"]*h
        assert mat["positive"](bare)
        if log_x == 0:
            graph = bare
        else:
            assert mat["positive"](A.inv()-bare)
            graph = (arb_mat([[1,0,0],[0,1,0],[0,0,1]])-bare*A).inv()*bare
        c, r = ports_at(-s)
        side = z+9*c+r
        leg = arb_mat([[-6*c], [2*r], [side]])*(1/(side*(b+s*d)+4))
        leg_x = (leg.transpose()*(family+A*graph*family))[0, 0]
        slope_b = (e*integral_slope(e/triplet**2)-s*integral_slope(s/triplet**2))
        slope_b *= 48/(16*arb.pi()**2*b*triplet**2)
        gradient = [-k*s*leg_x**2*(bj+s*dj)+graph[0,2]**2*slope_b*bj
                    +factor*(v.transpose()*graph*family)[0, 0]**2
                    for bj, dj, factor, v in zip(db, dd, factors, vectors, strict=True)]
        values.append(graph[2, 2])
        derivs.append(gradient)
    residual, jacobian = [], []
    for i, ref in enumerate(refs, 1):
        norm = (q(".5")*i).exp()*values[0]*ref
        residual.append(values[i]/norm-1)
        jacobian.append([(g-values[i]*o/values[0])/norm
                         for g, o in zip(derivs[i], derivs[0], strict=True)])
    return arb_mat([[v] for v in residual]), arb_mat(jacobian), b, d


def certify_root(data, center):
    center = list(map(q, center))
    residual, jac, _, _ = readouts(data, center)
    approximate = jac.inv()
    pre = arb_mat([[approximate[i,j].mid() for j in range(2)] for i in range(2)])
    correction = pre*residual
    center = [(v-correction[j,0]).mid() for j,v in enumerate(center)]
    residual, _, _, _ = readouts(data, center)
    radius = q("2e-12")
    box = [v+arb(0, radius.upper()) for v in center]
    _, jac_box, b, d = readouts(data, box)
    defect = eye-pre*jac_box
    norm = max(sum(abs(defect[i,j]) for j in range(2)).upper() for i in range(2))
    offset = max(abs(v).upper() for v in (pre*residual).entries())
    print("contraction", norm.str(6), "offset", offset.str(6), flush=True)
    assert not pre.det().contains(0) and norm < 1 and offset+radius*norm < radius
    assert all(v > 0 for v in box)
    return box, b, d


def count_modes(data, ports, strengths, b, d):
    energies, z, k, h, _, _, atoms, _ = data
    e, second = energies
    # At the cut, eliminate the old fields using a symmetric bordered matrix.
    cut = q(".11")
    assert e < atoms[0][0] < second < cut < atoms[1][0]
    (pc, vc), (pr, vr) = ([arb_poly(a) for a in pair] for pair in ports)
    ac = vc*z-pc*(9/(24*unit))
    ar = -pr*(1/(4*unit))
    old = (ac*vr+ar*vc)*arb_poly([b, -d*unit])+4*vc*vr
    assert old[0] == 0 and not old[1].contains(0)
    pos = arb_poly(old.coeffs()[1:])
    old_count = ordered["variations_in"](pos, q(0), cut/unit)
    assert old_count in (0, 1, 2)
    # Descartes' upper bound 2 alone would permit 0. Both cases below
    # need only an upper bound, since three distinct roots are constructed.
    assert not old(atoms[0][0]/unit).contains(0)
    bare = mat["massive_exchange"](-cut, z, b, d)*k+mat["hypercharge"]*h
    small = arb_mat([[bare[i,j] for j in range(2)] for i in range(2)])
    mixing = arb_mat([[scalar(v, small, w) for _, w in atoms] for _, v in atoms])
    for j, ((t, _), u) in enumerate(zip(atoms, strengths, strict=True)):
        mixing[j,j] += (t-cut)*(t+e)/((cut+e)*u)
    assert not mixing.det().contains(0)
    if mixing.det() < 0:
        negatives = 1
    else:
        assert not mixing[0,0].contains(0)
        negatives = 2 if mixing[0,0] < 0 else 0
    assert old_count+negatives == 3
    kinetic = arb_mat([[k/24, 0], [0, k/4]])+(h-k/24)*tensor(y)
    trace = sum(u*scalar(v, kinetic, v)/(t+e) for u, (t,v) in zip(strengths,atoms,strict=True))
    assert mat["positive"](arb_mat([[k/24,0,0],[0,k/4,0],[0,0,k/d]])+mat["hypercharge"]*(h-k/24))
    assert 0 < trace < q(".01")
    return trace


def verify():
    time = poles["study"]["nominal"][2]
    target = poles["study"]["comp"]["model"]["target"](11432, time)
    *_, colour, right = poles["neutral"]["graph_constants"](port_builder=poles["neutral"]["pole_polynomials"])
    ports = (colour, right)
    for i, logs in enumerate(poles["retune"]["log_sets"]):
        (L,z,rho,eta,c,masses),_,_ = poles["retune"]["predictions"](logs,time)
        b0, d0 = c/rho, eta/rho
        n1 = poles["reference_readings"](logs)[3]
        e = n1**2/L
        other = mat["geometry"]["second_neutral_mass"](L,z,b0,d0)**2/L
        h = mat["flow"]["running"]["hypercharge_window"](time+n1.log()).mid()**2
        k = 1/(poles["offset"]+target/(15*L))
        refs = retuned["readings"](e,z,b0,d0,masses[0],1/rho,h,q(0))
        atoms = [dark_atom(ports,h,k,masses[0]), (q(".15"),arb_mat([[0],[1]]))]
        data = (e,other), z,k,h,masses[0],b0,atoms,refs
        box, b, d = certify_root(data, centres[i])
        trace = count_modes(data, ports, box, b, d)
        shift = ((b/b0).sqrt()-1)*10**6
        assert q("33.6") < shift < q("34")
        assert 4*masses[0]**2*(b/b0) < atoms[0][0]
        print("case",i,"triplet mass shift (ppm)",shift.str(9),
              "UV trace",trace.str(9),"shape D/D0",(d/d0).str(10),flush=True)


if __name__ == "__main__":
    verify()
