"""Can a uniformly sampled nineteen-point Behrend shell supply the four color harmonics?"""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
from runpy import run_path
from collections import Counter
from flint import acb, arb, arb_mat, ctx

ROOT = Path(__file__).resolve().parent
rank = run_path(str(ROOT/"ch39-color-moment-rank.py"))
two, q = rank["two"], rank["q"]
SEEDS = {
    3200: (996648149856942738148966446382121113162578877162050, 996862799511113882178236260793316095606307695373385, 997018161993491440782141135440767370236937101631927, 997200965004013610428062595764621996374917244534713),
    11485: (636680776750357401674305837384122728266212091568484, 835997088820406194585609360449107788900049124606139, 872877876948854952112777203303467913078604741376133, 955827006139319765181522280724753335999330236972028),
    11486: (636608048975196417997977882176096128976134326298878, 835945690467361484657541684702101671349235535051108, 872864859901408879475903637937417062026083840678166, 955819211372846554072143921011036798978622318747617),
    11487: (636497006335029519452914144601789773392237242016728, 835867222943411891351454871850540806489972160427342, 872844983368314507119196389681068497757611235114002, 955807308611722463206704325455234963591590531095159),
    11488: (636449220958129254231991382109066438941808303527418, 835833458987364742467599499583479225747595861745094, 872836429172376832157076544201943465617792791555811, 955802185893028045980396323944993623488568593881580),
    11489: (636348273219603198644863274426038972059328869298220, 835762138224953101381707775843505434997264837260821, 872818356937198483533452500570334909829774129550859, 955791362921460503119876821640200039631715655224236),
    11490: (636300578845709355121440086557682295405880133730484, 835728444592464553490940262056166013934385996167253, 872809817811118672091678067797246029182244745880640, 955786248913460950886244106453172651700363023626006),
}


def fixed_weights(target, counts, center, refinements=0):
    n = rank["N"]
    x = q(rank["EPS"]/3)
    phase = acb((1-x*x)/(1+x*x),2*x/(1+x*x))
    size = len(counts)
    assert size == len(target) == len(center)
    ratios = [(1-(phase**j).real)/(1-phase.real) for j in range(1,size+1)]
    weights = [q(F(c,sum(counts))) for c in counts]
    def evaluate(point):
        roots = [v**(q(1)/n) for v in point]
        residual, jac = [], []
        for r, goal in zip(ratios,target):
            bases = [1-r*(1-t) for t in roots]
            weighted = [w*b**n for w,b in zip(weights,bases)]
            residual.append([sum(weighted)-goal])
            jac.append([v*r*t/(u*b) for v,t,u,b in zip(weighted,roots,point,bases)])
        return arb_mat(residual), arb_mat(jac)
    for _ in range(refinements):
        res, jac = evaluate(center)
        correction = jac.solve(res)
        center = [(v-d).mid() for v,d in zip(center,correction.entries())]
    radius = q("1e-33")
    box = [v+arb(0,radius) for v in center]
    assert 0 < box[0] and box[-1] < 1 and all(a < b for a,b in zip(box,box[1:]))
    assert 1-box[0]**(q(1)/n) < 1-phase.real  # each non-idle probability is below one
    res, mid = evaluate(center)
    _, whole = evaluate(box)
    inv = arb_mat(size,size,[v.mid() for v in mid.inv().entries()])
    assert not inv.det().contains(0)
    err = arb_mat([[int(i==j) for j in range(size)] for i in range(size)])-inv*whole
    defect = inv*res
    assert all(abs(defect[i,0])+radius*sum(abs(err[i,j]) for j in range(size)) < radius for i in range(size))


def atom_bounds(target, nodes):
    a,b = nodes
    # Previous three-band certificate ensures every allowed distribution charges all three intervals.
    rank["require_three"](target, nodes)
    assert q(".705") < a < q(".707") and q(".881") < b < q(".884")
    delta = q("2e-14")
    m = [q(1),target[0],*[v+arb(0,delta) for v in target[1:]]]
    tangent = lambda v: v**4-4*b**3*v+3*b**4
    average = m[2]-4*b**3*m[1]+3*b**4
    assert average/tangent(a) < q(F(3,19)) and average/tangent(q(0)) < q(F(1,19))
    # A nonnegative polynomial with two interior double roots, normalized at one.
    t, u = q(".66"), q(".875")
    rows = [rank["row"](v,d) for v,d in [(t,0),(t,1),(u,0),(u,1)]]
    coeff = arb_mat([r[:-1] for r in rows]).solve(arb_mat([[-r[-1]] for r in rows]))
    poly = [*coeff.entries(),q(1)]
    assert sum(poly)>0 and sum(c*v for c,v in zip(poly,m))/sum(poly) < q(F(1,19))


def exclude_four(counts, target):
    # The exactly known first moment eliminates the fourth node.
    weights = [q(F(c,19)) for c in counts]
    scales = (4,9,16)
    goals = [v+arb(0,q("2e-14")) for v in target[1:]]
    domain = q(0).union(q(1))
    stack = [[(q(0),q(".707")),(q(".705"),q(".884")),(q(".705"),q(1))]]
    count = 0
    while stack:
        bounds = stack.pop()
        count += 1
        xs = [a.union(b) for a,b in bounds]
        raw = (target[0]-sum(w*x for w,x in zip(weights,xs)))/weights[3]
        if not raw.overlaps(domain):
            continue
        last = raw.intersection(domain)
        nodes = xs+[last]
        if last < q(".881") or any(x>y for x,y in zip(nodes,nodes[1:])):
            continue
        if any(not sum(w*x**k for w,x in zip(weights,nodes)).overlaps(v) for k,v in zip(scales,goals)):
            continue
        centers = [x.mid() for x in xs]
        z = (target[0]-sum(w*x for w,x in zip(weights,centers)))/weights[3]
        mid = centers+[z]
        residual = arb_mat([[sum(w*x**k for w,x in zip(weights,mid))-v] for k,v in zip(scales,goals)])
        def jac(first, fourth):
            return arb_mat([[k*w*(x**(k-1)-fourth**(k-1)) for w,x in zip(weights,first)] for k in scales])
        try:
            inverse = jac(centers,z).inv()
        except ZeroDivisionError:
            inverse = None
        if inverse is not None:
            pre = arb_mat(3,3,[v.mid() for v in inverse.entries()])
            defect = arb_mat([[int(i==j) for j in range(3)] for i in range(3)])-pre*jac(xs,raw)
            offset = arb_mat([[x-c] for x,c in zip(xs,centers)])
            enclosure = arb_mat([[c] for c in centers])-pre*residual+defect*offset
            if any(not x.overlaps(y) for x,y in zip(xs,enclosure.entries())):
                continue
            xs = [x.intersection(y) for x,y in zip(xs,enclosure.entries())]
        i = max(range(3), key=lambda j: float(xs[j].rad()))
        assert xs[i].rad() > q("1e-8"), "Unresolved moment box"
        bounds = [(x.lower(),x.upper()) for x in xs]
        lo, hi = bounds[i]
        midpt = (lo+hi)/2
        for interval in ((lo,midpt),(midpt,hi)):
            child = bounds.copy()
            child[i] = interval
            stack.append(child)
    return count


def verify():
    ctx.prec = 360
    shell = [p for p in product(range(30), repeat=3) if sum(x*x for x in p)==867]
    assert Counter(map(sum,shell)) == {51:1,49:3,47:9,35:6}
    assert sum(a==b==c for a,b,c in shell)==1 and len(shell)==19
    mu = 1000*q(F(429,16)).exp()
    common = (two["geo"]["BOOK_PLANCK"]/mu)**2/(8*arb.pi())
    compositions = [(1,3*i,3*j,3*(6-i-j)) for i in range(1,5) for j in range(1,6-i)]
    for index, seed in SEEDS.items():
        power, degree = (q(144196081),1650) if index==3200 else (common,190)
        target = two["harmonics"](acb.zeta_zero(index).imag,power,degree,4)
        center = [q(F(x,10**51)) for x in seed]
        fixed_weights(target, (1,6,9,3) if index==3200 else (1,2,15,1), center)
        if index!=3200:
            nodes = [q(F(x,10**18)) for x in rank["SEEDS"][index][:2]]
            atom_bounds(target,nodes)
            total = sum(exclude_four(c,target) for c in compositions)
            print("Zero, excluded four-node boxes:", index, total, flush=True)
        else:
            print("The 3200th zero admits the four actual radial multiplicities.", flush=True)


if __name__ == "__main__":
    verify()
