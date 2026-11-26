"""两处族流读数对完整重多重态内质量分裂的必要条件。"""
from itertools import combinations_with_replacement, product
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat

dirac = run_path(str(Path(__file__).with_name("ch39-two-dirac-current-witness.py")))
shape, mat, poles, q = (dirac[key] for key in ("shape", "mat", "poles", "q"))
centres = ((".00030685", ".25294657"), (".00027904", ".25303578"), (".00015810", ".25342516"))


def minimal_multiplet():
    # Symmetric square of the antifundamental and spin 7/2 of SU(2)_R.
    d_weights = [a+b for a, b in combinations_with_replacement((-1, -1, -1, 3), 2)]
    weights = list(product(d_weights, range(-7, 8, 2), (-4,)))
    assert (6, 7, -4) in weights and len(weights) == 80
    charges = arb_mat(weights)
    assert charges.transpose()*charges == arb_mat([[576, 0, 0], [0, 1680, 0], [0, 0, 1280]])
    def dimension(a, b, c):
        return (a+1)*(b+1)*(c+1)*(a+b+2)*(b+c+2)*(a+b+c+3)//12
    assert min(dimension(3,0,0), dimension(0,3,0), dimension(0,0,3)) >= 10
    small = {labels for labels in product(range(3), repeat=3) if dimension(*labels) < 10}
    assert small == {(0,0,0), (1,0,0), (0,1,0), (0,0,1)}
    right_index = 10*sum(r*r for r in range(-7, 8, 2))//4
    b_min = q(4)*right_index/3-q(22)/3
    doubling_bound = (8*arb.pi()**2/(b_min*q(2).log())).sqrt()
    assert right_index == 420 and q(".453") < doubling_bound < q(".455")
    return charges.transpose()*charges


def verify():
    trace = minimal_multiplet()
    charge = dirac["mixed_charge"]
    companions = trace-charge*charge.transpose()
    assert mat["positive"](companions)
    time = poles["study"]["nominal"][2]
    target = poles["study"]["comp"]["model"]["target"](11430, time)
    for i, logs in enumerate(poles["retune"]["log_sets"]):
        (scale, z, rho, eta, c, masses), _, _ = poles["retune"]["predictions"](logs, time)
        neutral = poles["reference_readings"](logs)[3]
        e = neutral**2/scale
        h = mat["flow"]["running"]["hypercharge_window"](time+neutral.log()).mid()**2
        k = 1/(poles["offset"]+target/(15*scale))
        times = q("1/2"), q(1)
        graphs, targets, gaps, full_packets = [], [], [], []
        for log_x in times:
            p, a, _ = shape["probe"](log_x, e, z, c/rho, eta/rho, masses[0])
            g = ((p*k+mat["hypercharge"]*h).inv()-a).inv()
            assert mat["positive"](g)
            graphs.append(g)
            x, _, dest, gap, *_ = dirac["packet"](log_x, e, z, c/rho, eta/rho, masses[0], 1/rho, k, h)
            targets.append(dest)
            gaps.append(gap)
            vacuum = companions*dirac["dirac"](x, e/(400*masses[0])**2)[0]
            inverse = g.inv()-vacuum
            assert mat["positive"](inverse)
            full = inverse.inv()
            remaining = gap-(g*vacuum*full)[2, 2]
            vertex = full*charge
            full_packets.append((x, full[2, 2], dest, remaining, vertex[2, 0]**2,
                                 (charge.transpose()*vertex)[0, 0]))
        tolerance = q("2e-8")
        low, high = gaps[0]-tolerance*targets[0], gaps[1]+tolerance*targets[1]
        x1, x2 = (value.exp() for value in times)
        threshold = 4*masses[0]**2/e
        cap = (x1-1)/(x2-1)*(threshold+x2)/(threshold+x1)
        u, v = graphs
        # In a complete representation the charge trace is diagonal. For a
        # diagonal budget D, v[:,2]' D v[:,2] bounds each positive component.
        assert low > 0 and high > 0
        entry_ratio = max((u[j, 2]**2/v[j, 2]**2).upper() for j in range(3))
        trace_ratio = max((u[j, j]/v[j, 2]**2).upper() for j in range(3))
        mass_ratio = (low/(cap*high*(entry_ratio+low*trace_ratio))).sqrt()
        print("case", i, "mass ratio lower bound", mass_ratio.lower().str(10),
              "allowed relative amplitude budget at zero splitting",
              ((cap*high*entry_ratio/(1-cap*high*trace_ratio))/low).str(10), flush=True)
        assert mass_ratio > 1
        # At fixed r*u*(1-u) the ratio of the two logarithms increases
        # with that scalar, so bound it at r = e/M3² and u = 1/2.
        max_t = e/(4*masses[0]**2+e)
        fermion_cap = (1+(x1-1)*max_t).log()/(1+(x2-1)*max_t).log()
        ratios = [u[j, 2]/v[j, 2] for j in range(3)]
        hi = max(r.upper() for r in ratios)
        lo = min(r.lower() for r in ratios)
        radius = (hi-lo)/2
        splitting = q("28.8")
        projected = hi**2+hi*(splitting-1/splitting)*radius+splitting**2*radius**2
        budget = fermion_cap*high
        den = 1-splitting**2*budget*trace_ratio
        assert lo > 0 and den > 0 and budget*projected/den < low
        centre = [q(value) for value in centres[i]]
        jacobian = arb_mat([dirac["equation"](row, *centre)[1] for row in full_packets])
        fitted = jacobian.inv()
        preconditioner = arb_mat([[fitted[a,b].mid() for b in range(2)] for a in range(2)])
        boxes, contraction = dirac["certify"](full_packets, centre, preconditioner)
        mass1, mass2 = [(e/r).sqrt()/masses[0] for r in boxes]
        hierarchy = 400/mass2
        assert 1 < mass2 < mass1 < 400 and 353 < hierarchy < 355
        print("case", i, "28.8 splitting maximal fraction", (budget*projected/den/low).str(9),
              "companion construction mass ratios", mass1.str(8), mass2.str(8),
              "hierarchy", hierarchy.str(8), "contraction", contraction.str(8), flush=True)


if __name__ == "__main__":
    verify()
