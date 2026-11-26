"""Enclose each zero's input-volume-weighted score in the original low-input box."""
from dataclasses import dataclass
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, ctx

root = Path(__file__).parent
base = run_path(str(root / "ch39-low-input-volume-ranking.py"))
gauge, q, running, ratios = (base[k] for k in ("gauge", "q", "running", "ratios"))
loop = gauge["loop"]
anchor = q(gauge["fixed"]["L"])
radius = q("1/8")


def flow_jacobian(g, y):
    """Jacobian in the three inverse gauge and nine log-Yukawa coordinates."""
    matrix = arb_mat(12, 12)
    for i in range(3):
        for j in range(3):
            matrix[i, j] = 2*running["B"][i][j]*g[j]**2/loop**2
        for a in range(9):
            matrix[i, 3+a] = 2*running["D"][i][a//3]*y[a]/loop**2
    for a in range(9):
        family = a//3
        for j in range(3):
            matrix[3+a, j] = 2*running["C"][family][j]*g[j]**2/loop
        for b in range(9):
            mix = family < 2 and b == 3*(1-family)+a%3
            coefficient = 2*(3 if b < 6 else 1) + 3*(int(a == b)-int(mix))
            matrix[3+a, 3+b] = coefficient*y[b]/loop
    return matrix


@dataclass
class FlowBounds:
    faces: list
    top_sensitivity: list
    target_ratios: list
    transformed_rates: list
    rates: list
    determinant: arb
    residual_norms: list
    low_box: list
    initial_offsets: list


def variations(low, logs, *, refinements=0):
    """Integrate the full variational equation and enclose its first three rows on the energy band."""
    tangent = arb_mat([[int(i == j) for j in range(12)] for i in range(12)])
    block, rates = None, None
    elapsed = arb(0)

    def observe(step, g, y, traces):
        nonlocal tangent, block, rates, elapsed
        field = flow_jacobian(g, y)
        norm = max(sum(abs(field[i,j]) for j in range(12)).upper() for i in range(12))
        col = [max(abs(tangent[i,j]).upper() for i in range(12)) for j in range(12)]
        # For each column the integral equation bounds its in-step increment
        # by h ||J|| exp(h ||J||) times its initial maximum component.
        change = step*norm*(step*norm).exp()
        pad = arb_mat([[arb(0,change*v) for v in col]]*12)
        derivative = field*(tangent + pad)
        for _ in range(refinements):
            derivative = field*(tangent + arb(0).union(step)*derivative)
        if elapsed+step >= anchor-radius:
            slab = tangent + arb(0).union(step)*derivative
            rows = [[slab[i,j] for j in range(12)] for i in range(3)]
            flow_rate = [-2*(running["BETA"][i] + sum(b*u-d*t for b,d,u,t in
                         zip(running["B"][i], running["D"][i], g, traces))/loop)/loop
                         for i in range(3)]
            block = [[a.union(b) for a,b in zip(old,new)] for old,new in zip(block,rows)] if block else rows
            rates = [a.union(b) for a,b in zip(rates,flow_rate)] if rates else flow_rate
        tangent += step*derivative
        elapsed += step

    gauge["evolve"](low, logs, anchor+radius, 4096, diagonal_quarks=True, observe_flavors=observe)
    return arb_mat(block), rates


def prepare(gauge_lo, gauge_hi, logs):
    low = [a.union(b) for a,b in zip(gauge_lo,gauge_hi)]
    x_ref = [(a+b)/2 for a,b in zip(gauge_lo,gauge_hi)]
    y_ref = [y.mid() for y in logs]
    jac, rates = variations(low, logs)
    m = arb_mat([[jac[i,j] for j in range(3)] for i in range(3)])
    pre = arb_mat([[m[i,j].mid() for j in range(3)] for i in range(3)]).inv()
    assert m.det() > 0
    end, _, _ = gauge["evolve"](x_ref, y_ref, anchor, 65536, diagonal_quarks=True)

    def column(values):
        return arb_mat([[v] for v in values])

    normalized = pre*column(ratios)
    slope = pre*column(rates)
    center = pre*column(end)
    remainder = pre*jac
    for i in range(3):
        remainder[i,i] -= 1
    inputs = [v-x for v,x in zip(low,x_ref)] + [y-ref for y,ref in zip(logs,y_ref)]
    top = 5 # third up-type Yukawa in the full twelve-coordinate input
    residual = [sum(remainder[i,j]*inputs[j] for j in range(12) if j != top) for i in range(3)]
    faces = [(center[i,0]+residual[i]+(gauge_hi[i]-x_ref[i]),
              center[i,0]+residual[i]+(gauge_lo[i]-x_ref[i])) for i in range(3)]
    return FlowBounds(faces, [remainder[i,top] for i in range(3)], list(normalized),
                      list(slope), rates, m.det(),
                      [sum(abs(pre[i,j])*ratios[j] for j in range(3)) for i in range(3)],
                      low, [x_ref[i]-center[i,0]-residual[i] for i in range(3)])


def widths(target, data, top, tolerance=0, gauge_prior=None):
    # High and low faces of each input coordinate bound its first and last crossing.
    both = [[], []]
    for i in range(3):
        norm = data.target_ratios[i]
        assert norm > 0
        # Uniformly enclose the translated residual target for every point
        # in the three-dimensional tolerance cube, including all entry faces.
        shift = arb(0, tolerance*data.residual_norms[i])
        assert data.rates[i]/ratios[i] + 2*target*(-2*radius).exp() > 0
        for side, value in enumerate(data.faces[i]):
            start = (value + data.top_sensitivity[i]*top + shift)/norm
            rate = data.transformed_rates[i]/norm
            delta = (-radius).union(radius)
            assert target*(2*radius).exp() > start-radius*rate
            assert target*(-2*radius).exp() < start+radius*rate
            assert rate + 2*target*(-2*radius).exp() > 0
            for _ in range(6):
                delta = (target/(start+rate*delta)).log()/2
            both[side].append(delta)
    early, late = both
    lo = max(arb(0), min(x.lower() for x in late)-max(x.upper() for x in early)).lower()
    earliest, latest = max(x.lower() for x in early), min(x.upper() for x in late)
    if latest <= earliest:
        return [(arb(0), arb(0))]*2
    def densities(time):
        common = target*(-2*time).exp()
        speeds = [2*common + v/r for v,r in zip(data.rates,ratios)]
        assert all(speed > 0 for speed in speeds)
        if gauge_prior is None:
            return sum(speeds)/data.determinant, 1/data.determinant
        factor = 1/data.determinant
        initial = []
        for i, box in enumerate(data.low_box):
            x = (data.initial_offsets[i] + data.target_ratios[i]*common
                 - data.transformed_rates[i]*time - data.top_sensitivity[i]*top
                 + arb(0,tolerance*data.residual_norms[i]))
            if not box.overlaps(x):
                return None
            initial.append(box.intersection(x))
        factor *= gauge_prior(initial)
        assert factor > 0
        return (sum(speeds)*factor, factor)

    def integral(left, right, pieces, lower):
        if right <= left:
            return [arb(0),arb(0)]
        length = (right-left)/pieces
        total = [arb(0),arb(0)]
        for j in range(pieces):
            slab = (left+j*length).union(left+(j+1)*length)
            vals = densities(slab)
            if vals is None:
                assert not lower
                continue
            for k,value in enumerate(vals):
                total[k] += length*(value.lower() if lower else value.upper())
        return [(v.lower() if lower else v.upper()) for v in total]

    if gauge_prior is None:
        vals = densities(earliest.union(latest))
        return [((v*lo).lower(), (v*(latest-earliest)).upper()) for v in vals]
    minimum = integral(max(x.upper() for x in early), min(x.lower() for x in late), 16, True)
    maximum = integral(earliest, latest, 16, False)
    return list(zip(minimum,maximum))


def input_bounds():
    gauge_lo = [1/(c*q("1.001"))**2 for c in gauge["centers"]]
    gauge_hi = [1/(c*q(".999"))**2 for c in gauge["centers"]]
    logs = [2*(c*q(".99").union(q("1.01"))).log() for c in gauge["yukawa_centers"]]
    return gauge_lo, gauge_hi, logs


def zero_targets(tolerance):
    fixed = gauge["fixed"]
    common = (fixed["geo"]["BOOK_PLANCK"]/1000)**2 * 3*arb(5).sqrt()/(4*arb.pi())*(-2*anchor).exp()
    offset = q(fixed["OFFSET"])
    # Regenerate the necessary band over 23 <= L <= 33 with its independent enclosure.
    window = run_path(str(root / "ch39-floating-window-refinement.py"))
    outer, _ = window["bands"](window["scale_bounds"](tolerance), tolerance)
    first, last = window["base"]["certified_labels"](*outer)
    if not tolerance:
        assert (first,last) == (10815,12092)
    print("normalized residual half-width", tolerance, "outer labels", first, last)
    zeros = acb.zeta_zeros(first, last-first+1)
    assert len(zeros) == last-first+1 and all(z.real == q("1/2") for z in zeros)
    targets = {n:(common/(z.imag**2+offset),base["spectral"]["weight"](z.imag))
               for n,z in enumerate(zeros, first)}
    return targets, zeros[11737-first].imag


def verify(tolerance=0):
    ctx.prec = 112
    gauge_lo, gauge_hi, logs = input_bounds()
    data = prepare(gauge_lo, gauge_hi, logs)
    targets, height = zero_targets(tolerance)
    y0 = 2*(gauge["yukawa_centers"][2]*q(".99")).log()
    y1 = 2*(gauge["yukawa_centers"][2]*q("1.01")).log()
    measures = ("log Yukawa", "linear Yukawa")
    protocols = ("adjustable scale", "uniform scale")
    scores = {(protocol, measure):{} for protocol in protocols for measure in measures}
    candidates = list(targets)
    for pieces in (1, 32, 256):
        edges = [y0+(y1-y0)*j/pieces for j in range(pieces+1)]
        bins = [a.union(b)-logs[2].mid() for a,b in zip(edges,edges[1:])]
        linear = [(b/2).exp() for b in edges]
        masses = [[arb(1)/pieces]*pieces,
                  [(b-a)/(linear[-1]-linear[0]) for a,b in zip(linear,linear[1:])]]
        for n in candidates:
            target, weight = targets[n]
            intervals = [widths(target,data,y,tolerance) for y in bins]
            for mode, protocol in enumerate(protocols):
                for measure, mass in zip(measures,masses):
                    lo = weight*sum(p*val[mode][0] for p,val in zip(mass,intervals))
                    hi = weight*sum(p*val[mode][1] for p,val in zip(mass,intervals))
                    scores[protocol, measure][n] = lo.lower(), hi.upper()
        survivors = set()
        for rankings in scores.values():
            best = max(v[0] for v in rankings.values())
            survivors.update(n for n,v in rankings.items() if not v[1] < best)
        candidates = sorted(survivors)
        print(pieces, "top-Yukawa bins, possible winners", candidates)
    assert candidates == [11737]
    for (protocol,measure), rankings in scores.items():
        lo, hi = rankings[11737]
        rival = max(v[1] for n,v in rankings.items() if n != 11737)
        margin = lo/rival - 1
        assert margin > 0
        if not tolerance:
            assert margin > q(".0010946" if protocol == "adjustable scale" else ".0021401")
        elif tolerance == arb(2)**-21:
            assert margin > q(".0005196" if protocol == "adjustable scale" else ".0015648")
        if not tolerance and (protocol,measure) == ("adjustable scale","log Yukawa"):
            assert lo > q("3.07870e-10") and rival < q("3.07535e-10")
        print(protocol, measure, "11737 score", lo.str(9), hi.str(9),
              "largest rival upper bound", rival.str(9), "relative advantage >", margin.lower().str(9))
    curvature = 1 + (1+(height/312)**2).sqrt()
    assert q("11346.10207491353") < height < q("11346.10207491355")
    assert q("37.37945839566") < curvature < q("37.37945839567")


if __name__ == "__main__":
    verify()
