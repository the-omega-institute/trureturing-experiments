"""Exhaust all Gaussian-selected cyclic orbit rays for zeta zeros 10 and 527."""
from itertools import combinations
from pathlib import Path
import runpy
from flint import acb, acb_mat, acb_poly, arb, arb_poly

fit = runpy.run_path(str(Path(__file__).with_name("ch39-hierarchy-in-zero-flux.py")))
src, block, wall, interior = (fit[key] for key in ("src", "block", "W", "J"))
nodes = [arb(j)/3 for j in range(4)]
identity = acb_mat([[int(i == j) for j in range(3)] for i in range(3)])


def adjoint(mat):
    return mat.conjugate().transpose()


def representative(face):
    sites = fit["s"]["sites"]
    cycles = []
    for _ in range(3):
        cycles.extend(face[j:] + face[:j] for j in range(len(face)))
        face = tuple(sites.index(sites[v][1:] + sites[v][:1]) for v in face)
    return min(cycles)


def wall_samples(m, orbit):
    c = block(m, interior, interior).inv()
    projection = block(m, wall, interior)*c
    k = block(m, wall, wall) - projection*block(m, interior, wall)
    h = acb_mat([[v[j,0] for v in orbit] for j in range(19)])
    hi = block(h, interior, range(3))
    v = block(h, wall, range(3)) - projection*hi
    g = adjoint(hi)*c*hi
    assert all(block(g, range(n), range(n)).det().real > 0 for n in range(1,4))
    # x=s/(4+s). Clearing det(L)>0 leaves all bond phases and singular ratios invariant.
    samples = []
    for x in nodes:
        L = (1-x)*identity+4*x*g
        samples.append(-L.det().real*(k+4*x*v*L.solve(adjoint(v))))
    return h, samples


def polynomials(number):
    height = acb.zeta_zero(number).imag
    m, orbits = fit["model"](*(fit["s"][key]*height for key in ("a","b")))
    assert all(block(m, range(n), range(n)).det().real > 0 for n in range(1,20))
    scores, loops, spectra = [], [], []
    for orbit in orbits:
        h, samples = wall_samples(m, orbit)
        resistance = adjoint(h)*m.solve(h)
        det_values = [((1-x)*identity+4*x*resistance).det().real for x in nodes]
        scores.append(arb_poly.interpolate(nodes, det_values))
        def bond(i,j):
            return acb_poly.interpolate(nodes,[mat[wall.index(i),wall.index(j)] for mat in samples])
        current = {}
        for f in src["all_faces"]():
            poly = acb_poly([1])
            for i,j in zip(f, f[1:]+f[:1]):
                poly *= bond(i,j)
            current[f] = poly
        loops.append(current)
        symbols = [fit["symbols"](block(mat,fit["ii"],fit["oo"])) for mat in samples]
        spectra.append([acb_poly.interpolate(nodes,[row[j] for row in symbols]) for j in range(3)])
    return scores, loops, spectra


def isolate(poly, real_part=None):
    """Isolate all simple real roots on [0,1], optionally restricted to Re<0."""
    result = []
    derivative = poly.derivative()
    todo = [(arb(0),arb(1))]
    while todo:
        a,b = todo.pop()
        interval = a.union(b)
        if not poly(interval).contains(0) or (real_part is not None and real_part(interval)>0):
            continue
        if not derivative(interval).contains(0):
            ends = poly(a)*poly(b)
            if ends>0:
                continue
            if ends<0 and (real_part is None or real_part(interval)<0):
                while b-a>arb("1e-17"):
                    mid=(a+b)/2
                    sides=poly(a)*poly(mid)
                    if sides<0:
                        b=mid
                    else:
                        assert sides>0
                        a=mid
                result.append(a.union(b))
                continue
        assert b-a>arb("1e-18")
        mid=(a+b)/2
        todo.extend([(a,mid),(mid,b)])
    return result


def cuts(loops):
    events = []
    seen = set()
    for f, poly in loops.items():
        key = representative(f)
        if key in seen:
            continue
        seen.add(key)
        real, imag = (arb_poly([getattr(z,part) for z in poly.coeffs()]) for part in ("real","imag"))
        events.extend((root, f) for root in isolate(imag,real))
    return events


def divide_axis(events):
    events.sort(key=lambda item:float(item[0].mid()))
    ends = [arb(0)] + [root for root,_ in events] + [arb(1)]
    assert all(a<b for a,b in zip(ends,ends[1:]))
    return ends


def ratios_above(symbols, interval, ratio):
    vals = [abs(poly(acb(interval))) for poly in symbols]
    upper = max(v.upper() for v in vals)
    return sum(v.lower()>ratio*upper for v in vals)>=2


def check(number, wanted, floor):
    scores, loops, spectra = polynomials(number)
    ties=[]
    for i,j in combinations(range(4),2):
        delta=scores[i]-scores[j]
        assert delta[0].contains(0)  # all relative determinants equal one at s=0
        ties.extend((root,(i,j)) for root in isolate(delta.right_shift(1)))
    ends=divide_axis(ties)
    def winner(x):
        vals=[score(x) for score in scores]
        idx=min(range(4),key=lambda i:float(vals[i].mid()))
        assert all(vals[idx]<v for i,v in enumerate(vals) if i!=idx)
        return idx
    winners=[winner((a.upper()+b.lower())/2) for a,b in zip(ends,ends[1:])]
    # Check every crossing itself including ties that do not change the minimizer.
    for (root,pair),a,b in zip(ties,winners,winners[1:]):
        if a==b:
            assert all(scores[a](root)<score(root) for i,score in enumerate(scores) if i!=a)
        else:
            assert {a,b}==set(pair)
            assert all(scores[a](root)<score(root) for i,score in enumerate(scores) if i not in pair)
    windows=[]
    for orbit in set(winners):
        events=cuts(loops[orbit])
        axis=divide_axis(events)
        polys={f:poly.coeffs() for f,poly in loops[orbit].items()}
        for left,right in zip(axis,axis[1:]):
            deg=src["degrees"](polys,(left.upper()+right.lower())/2)
            if deg!=(3,[3]):
                continue
            for i,winning in enumerate(winners):
                if winning!=orbit or right<ends[i] or ends[i+1]<left:
                    continue
                # These examples have flux crossings strictly inside a single action chamber.
                assert ends[i]<left and right<ends[i+1]
                assert ratios_above(spectra[orbit],left.union(right),arb(floor))
                windows.append((orbit, left,right))
    assert len(windows)==1 and windows[0][0]==wanted
    orbit,left,right=windows[0]
    print("zero, orbit, bare stiffness endpoints, lower mass ratio",number,orbit,
          *( (4*x/(1-x)).str(15) for x in (left,right)),floor)
    return windows[0]


if __name__=="__main__":
    check(10,3,".22")
    check(527,2,".20")
