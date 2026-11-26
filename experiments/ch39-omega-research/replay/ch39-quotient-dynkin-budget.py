"""公共中心商的矢量型物质在两项族流读数下的最小颜色与右弱指数之和。"""
from collections import Counter
from fractions import Fraction
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import arb_mat

model = run_path(str(Path(__file__).with_name("ch39-two-dirac-current-witness.py")))
shape, mat, poles, q = (model[k] for k in ("shape", "mat", "poles", "q"))
mass46 = (".28549878",".28532298",".28456476")
mass60 = ((".23280505",".29949121"),(".23239576",".29973955"),(".23061972",".30081799"))


def colour_weights(a, b, c):
    # SU(4) -> U(3) interlacing, including each SU(3) multiplicity.
    first, second = a+b+c, b+c
    total = first+second+c
    weights = Counter()
    for i,j,k in product(range(second, first+1), range(c, second+1), range(c+1)):
        dim = (i-j+1)*(j-k+1)*(i-k+2)//2
        weights[4*(i+j+k)-3*total] += dim
    return weights


def content(a, b, c, p, left_dim):
    assert (a+2*b+3*c+p+left_dim-1)%2 == 0
    weights = {(d,r): n*left_dim for d,n in colour_weights(a,b,c).items()
               for r in range(-p,p+1,2)}
    colour = Fraction(sum(d*d*n for (d,_),n in weights.items()), 24)
    right = Fraction(sum(r*r*n for (_,r),n in weights.items()), 4)
    return weights, colour, right


def catalogue(twice_budget):
    # The Weyl dimension and quadratic Casimir are coordinatewise increasing.
    for labels in ((6,0,0), (0,6,0), (0,0,6)):
        weights = colour_weights(*labels)
        assert sum(d*d*n for d,n in weights.items()) > 12*twice_budget
    assert 12*13*14 > 6*twice_budget
    entries = []
    for a,b,c,p in product(range(6), range(6), range(6), range(12)):
        if not (a or b or c or p) or c>a:
            continue
        weights, colour, right = content(a,b,c,p,1+(a+2*b+3*c+p)%2)
        cost = 2*(colour+right)
        if cost <= twice_budget:
            assert cost.denominator == 1
            entries.append((int(cost),weights))
    return entries


def maximal_pairs(pairs):
    # Keep only coordinatewise maximal pairs of integer upper bounds.
    highest = -1
    frontier = []
    for a,b in sorted(set(pairs),reverse=True):
        if b > highest:
            frontier.append((a,b))
            highest = b
    return frontier


def exclude(data, graphs, cap_mass, twice_budget, tolerance, angles, hypercharge=None):
    u, v = graphs
    low = data[0][3]-tolerance*data[0][2]
    high = data[1][3]+tolerance*data[1][2]
    x1, x2 = [row[0] for row in data]
    t = cap_mass/(4+cap_mass)
    cap = (1+(x1-1)*t).log()/(1+(x2-1)*t).log()
    fmax = model["dirac"](x1,cap_mass)[0]
    axes = arb_mat([[1,0],[0,1],[-v[0,2]/v[2,2],-v[1,2]/v[2,2]]])
    gram = axes.transpose()*u*axes
    budget = q(twice_budget)/2
    lam = ((budget*fmax*max((24*gram[0,0]).upper(),(4*gram[1,1]).upper())).sqrt()
           +(u[2,2]*cap*high).sqrt()/v[2,2])**2
    assert 0 < lam < 1 and low > 0
    entries = catalogue(twice_budget)
    weights = set().union(*(row.keys() for _,row in entries))
    scale = 10000
    xmax = 48
    scores = {}
    if hypercharge:
        e, triplet, penalty, log_budget = hypercharge
        fraction_at_q = model["dirac"](x1,e)[0]/fmax
        max_log = -triplet.log()
        assert 0 < fraction_at_q < 1 and penalty*log_budget < 1
    for angle in angles:
        ratio = angle/cap
        threshold = (1-lam)*low-angle*high
        assert threshold > 0
        quadratic = u[2,2]**2-ratio*v[2,2]**2
        assert quadratic < 0
        coeff = fmax*scale/threshold
        for d,r in weights:
            first, second = (d*g[0,2]+r*g[1,2] for g in graphs)
            lin = 2*(u[2,2]*first-ratio*v[2,2]*second)
            const = first**2-ratio*second**2
            assert quadratic*xmax*xmax+abs(lin)*xmax+const < 0
            assert 2*quadratic*xmax+abs(lin) < 0
            for charge in range(-xmax,xmax+1):
                z = (quadratic*charge+lin)*charge+const
                score = coeff*z
                if hypercharge:
                    charge_y = q(d)/6+q(r)/2
                    score = max((fraction_at_q*score).upper(),
                                (score-scale*penalty*charge_y**2*max_log).upper())/(1-penalty*log_budget)
                upper = max(0,int(score.upper().ceil().unique_fmpz()))
                scores.setdefault((d,r,charge),[]).append(upper)
    # Every irrep has one shared X. Round upward before grouping its weights.
    by_cost = {}
    for cost, mult in entries:
        options = by_cost.setdefault(cost,set())
        for charge in range(-xmax,xmax+1):
            val = tuple(sum(n*scores[d,r,charge][i] for (d,r),n in mult.items())
                        for i in range(2))
            options.add(val)
    items = [(c,a,b) for c,opts in by_cost.items() for a,b in maximal_pairs(opts)]
    fronts = [[(0,0)]]
    for capacity in range(1,twice_budget+1):
        new = list(fronts[-1])
        for cost,a,b in items:
            if cost <= capacity:
                new.extend((min(scale,x+a),min(scale,y+b)) for x,y in fronts[capacity-cost])
        fronts.append(maximal_pairs(new))
    peak = max(min(x,y) for x,y in fronts[-1])
    print("excluded half-index budget",twice_budget,"joint upper score",peak,"/",scale,flush=True)
    assert peak < scale
    return peak


def current_fields(total):
    # (a,b,c,p,left dimension,integer X,number of Dirac copies).
    fields = ([(2,0,0,2,1,-1,1), (0,1,0,2,1,-1,1), (1,0,0,0,2,-1,2)]
              if total==46 else [(2,0,0,2,1,-1,2), (0,0,0,2,1,-2,1)])
    spectrum = Counter()
    cs = rs = 0
    for a,b,c,p,left,x,n in fields:
        weights, colour, right = content(a,b,c,p,left)
        cs += n*colour
        rs += n*right
        for (d,r),num in weights.items():
            spectrum[d,r,x] += num*n
    assert cs+rs == total
    return spectrum


def fit_spectra(case, total, graphs, packets, cap_mass):
    spectrum = current_fields(total)
    qa = (-3,0,-1) if total==46 else (0,2,-2)
    qb = (2,2,-1)
    na, nb = spectrum.pop(qa), spectrum.pop(qb)
    a,b = [arb_mat([[entry] for entry in v]) for v in (qa,qb)]
    companions = arb_mat(3,3)
    for charge,n in spectrum.items():
        v = arb_mat([[w] for w in charge])
        companions += v*v.transpose()*n
    basis = arb_mat([[a[i,0], b[i,0]] for i in range(3)])
    fixed = []
    for g, (x,base,target,gap,*_) in zip(graphs,packets,strict=True):
        vacuum = companions*model["dirac"](x,cap_mass/10**6)[0]
        gv = g*vacuum
        assert sum(gv[i,i] for i in range(3)) < 1
        eye = arb_mat([[1,0,0],[0,1,0],[0,0,1]])
        full = (eye-gv).inv()*g
        remaining = gap-(gv*full)[2,2]
        projected = basis.transpose()*full
        vertex = arb_mat([[projected[i,2]] for i in range(2)])
        fixed.append((x,target,remaining,vertex,projected*basis))

    def equation(row, ra, rb):
        x,target,gap,vertex,gram = row
        fa, da = model["dirac"](x,ra)
        fb, db = model["dirac"](x,rb)
        factors = arb_mat([[na*fa,0],[0,nb*fb]])
        schur = arb_mat([[1,0],[0,1]])-gram*factors
        assert schur[0,0]>0 and schur.det()>0
        response = schur.inv()*vertex
        gain = (vertex.transpose()*factors*response)[0,0]
        grad = [response[0,0]**2*na*da/target, response[1,0]**2*nb*db/target]
        return (gain-gap)/target, grad

    if total==46:
        ra = q(mass46[case])
        assert 0 < ra < cap_mass
        errors = [equation(row,ra,cap_mass)[0]*10**9 for row in fixed]
        print("index 46 residuals ppb", *(v.str(8) for v in errors),flush=True)
        assert all(abs(v) < q("18.1") for v in errors)
        return
    centre = [q(v) for v in mass60[case]]
    values = [equation(row,*centre) for row in fixed]
    midpoint = arb_mat([g for _,g in values]).inv()
    pre = arb_mat([[midpoint[i,j].mid() for j in range(2)] for i in range(2)])
    rad = q(".000001")
    box = [c+model["arb"](0,rad.upper()) for c in centre]
    assert all(0 < v < cap_mass for v in box) and not pre.det().contains(0)
    jac = arb_mat([equation(row,*box)[1] for row in fixed])
    defect = arb_mat([[1,0],[0,1]])-pre*jac
    norm = max(sum(abs(defect[i,j]) for j in range(2)).upper() for i in range(2))
    shift = pre*arb_mat([[value] for value,_ in values])
    assert max(abs(s).upper() for s in shift.entries())+rad*norm < rad
    print("index 60 exact mass ratios",*((cap_mass/v).sqrt().str(8) for v in box),
          "contraction",norm.str(6),flush=True)


def verify():
    time = poles["study"]["nominal"][2]
    goal = poles["study"]["comp"]["model"]["target"](11430,time)
    for case, logs in enumerate(poles["retune"]["log_sets"]):
        (scale,z,rho,eta,c,masses),_,_ = poles["retune"]["predictions"](logs,time)
        neutral = poles["reference_readings"](logs)[3]
        e = neutral**2/scale
        h = mat["flow"]["running"]["hypercharge_window"](time+neutral.log()).mid()**2
        k = 1/(poles["offset"]+goal/(15*scale))
        graphs, packets = [], []
        for log_x in (q(".5"),q(1)):
            p, a, _ = shape["probe"](log_x,e,z,c/rho,eta/rho,masses[0])
            graph = ((p*k+mat["hypercharge"]*h).inv()-a).inv()
            assert mat["positive"](graph)
            graphs.append(graph)
            packets.append(model["packet"](log_x,e,z,c/rho,eta/rho,masses[0],1/rho,k,h))
        cap_mass = e/masses[0]**2
        exclude(packets,graphs,cap_mass,91,q("2e-8"),(q(".282"),q(".262")))
        exclude(packets,graphs,cap_mass,92,q("2e-8"),(q(".278"),q(".282")),
                hypercharge=(e,masses[0],q(".0213"),q(36)))
        for total in (46,60):
            fit_spectra(case,total,graphs,packets,cap_mass)


if __name__ == "__main__":
    verify()
