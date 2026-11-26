"""Phase, positive zeta height and original shell mass from a local return amplitude."""
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx
from sympy import fraction

port = run_path(str(Path(__file__).with_name("ch39-cayley-port-readout.py")))
shell = port["shell"]
mass = shell["MASS"]
ball = port["ball"]


def triangle_terms(neighbors, root):
    constant, linear = F(0), F(0)
    counts = Counter()
    for i, j in combinations(neighbors[root], 2):
        if j not in neighbors[i]:
            continue
        a, h = neighbors[root][i]
        b, k = neighbors[i][j]
        c, l = neighbors[j][root]
        winding = abs(h+k+l)
        counts[winding] += 1
        assert winding in (0, 1)
        if winding:
            linear += 2*a*b*c
        else:
            constant += 2*a*b*c
    return constant, linear, counts


def verify():
    ctx.prec = 192
    nodes, edges = shell["geometry"]()
    neighbors = [{} for node in nodes]
    for i, j, weight, charge in edges:
        neighbors[i][j] = weight, charge
        neighbors[j][i] = weight, -charge
    coeffs = [triangle_terms(neighbors, i) for i in range(len(nodes))]
    assert sum(b > 0 for _, b, _ in coeffs[:19]) == 12
    assert sum(b > 0 for _, b, _ in coeffs[19:]) == 6
    index = nodes.index((33, 33, 75))
    constant, slope, counts = coeffs[index]
    assert counts == {0: 12, 1: 4} and slope == F(3517425, 42342346661888)
    assert slope == max(b for _, b, _ in coeffs)
    cuts = [ball(constant+slope*F(q*q-1, q*q+1)) for q in port["cut_denominators"]]
    numerator, denominator = fraction(mass["response"](3))
    margins, masses = [], []
    for n in range(11485, 11491):
        rho = acb.zeta_zero(n)
        matrix = shell["matrices"](nodes, edges, (rho-1)/rho)
        returned = (matrix**3)[index, index]
        assert returned.imag.contains(0)
        reading = returned.real
        x = (reading-ball(constant))/ball(slope)
        height_squared = (1+x)/(4*(1-x))
        assert height_squared.overlaps(rho.imag**2) and height_squared.rad() < arb("1e-28")
        assert all(abs(reading-cut) > arb("1.19e-20") for cut in cuts)
        assert 11485 + sum(reading > cut for cut in cuts) == n
        margins.extend(abs(reading-cut) for cut in cuts)
        predicted = mass["evaluated"](numerator,x)/mass["evaluated"](denominator,x)
        old = acb_mat([[matrix[i,j] for j in range(19)] for i in range(19)])
        target, _ = shell["trace_values"](old)
        assert predicted.overlaps(target) and predicted.rad() < arb("1e-40")
        masses.append(predicted)
        print(n, "local return", reading.str(23), "original shell mass", predicted.str(25))
    assert all(b-a > arb("7e-24") for a, b in zip(masses, masses[1:]))
    print("minimum return margin", min(margins).str(14), "sensitive sites", 12, "+", 6)


if __name__ == "__main__":
    verify()
