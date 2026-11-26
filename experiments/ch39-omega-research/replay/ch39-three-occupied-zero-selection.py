"""Individually select six zeta heights by three occupied levels on a fixed positive graph."""
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx
from sympy import expand, symbols

shell = run_path(str(Path(__file__).with_name("ch39-refined-mass-obstruction.py")))
ball = shell["ball"]
weights = (28929, 28931, 28934, 28935, 28937, 28939)


def triangle_coordinate(x):
    # On [.99,1], 4c^3-3c increases strictly. Endpoints stay exact dyadics.
    lo, hi = F(63,64), F(1)
    assert ball(4*lo**3-3*lo) < x < 1
    for _ in range(165):
        c = (lo+hi)/2
        test = ball(4*c**3-3*c)
        if test < x:
            lo = c
        else:
            assert test > x
            hi = c
    return ball(lo).union(ball(hi))


def score(coordinate, weight):
    return (2*weight+1)*coordinate + (3*(1-coordinate**2)).sqrt()


def verify():
    ctx.prec = 224
    c, y = symbols("c y")
    assert expand((y-2*c)*((y+c)**2-3*(1-c**2))) == y**3-3*y-2*(4*c**3-3*c)
    lower = F(63,64)
    assert 3*(1-lower**2) < F(1,4) and min(weights)*(lower-F(1,2)) > 4
    nodes, edges = shell["geometry"]()
    triangle = [(21,51,69), (33,33,75), (39,39,69)]
    rotated = [(y,z,x) for x,y,z in triangle]
    marks = [[nodes.index(p) for p in group] for group in (triangle, rotated)]
    assert set(marks[0]).isdisjoint(marks[1])
    edge_charge = {(i,j): k for i,j,_,k in edges}
    edge_charge.update({(j,i): -k for i,j,_,k in edges})
    for group in marks:
        k = sum(edge_charge[i,j] for i,j in zip(group, group[1:]+group[:1]))
        assert abs(k) == 1
    eps = F(1, 10**20)
    degree = [0]*len(nodes)
    for i,j,_,_ in edges:
        if not any(i in g and j in g for g in marks):
            degree[i] += 1
            degree[j] += 1
    # Comparing sums of three eigenvalues at two heights costs at most 6||E||.
    disparity = ball(6*max(degree)*eps)
    assert disparity < arb("1e-17")
    coords = []
    for n in range(11484, 11492):
        rho = acb.zeta_zero(n)
        coords.append(triangle_coordinate(((rho-1)/rho).real))
    assert all(a < b for a,b in zip(coords,coords[1:]))
    for i, weight in enumerate(weights, 1):
        # In the actual graph the six marked edges and all 76 diagonal
        # entries stay fixed. The uniform Weyl bound holds at every zero.
        peak = score(coords[i], weight)
        margin = min(peak-score(coords[i-1], weight), peak-score(coords[i+1], weight))
        assert margin-disparity > arb("2.53e-14")
        print(11484+i, "integer triangle weight", weight, "global margin >", (margin-disparity).lower().str(10))


if __name__ == "__main__":
    verify()
