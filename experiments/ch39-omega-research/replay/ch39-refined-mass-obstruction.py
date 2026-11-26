"""Trace-convexity obstruction for a shared routing and mass graph."""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx
from sympy import cancel, diff, fraction, symbols

ROOT = Path(__file__).parent
LOCAL = run_path(str(ROOT/"ch39-local-shell-return.py"))
MASS = run_path(str(ROOT/"ch39-weighted-arc-zero-phase.py"))


def ball(q):
    value = F(q)
    return arb(value.numerator)/value.denominator


def geometry():
    pole = (51,51,51)
    sites = LOCAL["shell_grid"](3)|{pole}
    old = sorted({LOCAL["scaled"](p,3) for p in LOCAL["shell_grid"](1)}|{pole})
    nodes = old+sorted(sites-set(old))
    index = {p:i for i,p in enumerate(nodes)}
    squared = LOCAL["squared"]
    all_links = LOCAL["graph"](3,56)
    charge = {(p,q):w for p,row in all_links.items() for q,w in row}
    annulus = {p:[(q,w) for q,w in row if squared(p,q)>=108] for p,row in all_links.items()}
    endpoints = set(old)-{pole}
    minimum, paths = LOCAL["routes"](annulus, endpoints, True, 9)
    assert set(paths)==endpoints and max(minimum.values())==9
    row_sum = [F(0)]*len(nodes)
    edges, crossing, bridge_norm = [], 0, F(0)
    distances = [squared(p,q) for p,q in combinations(nodes,2)]
    assert min(distances)==12 and F(108,min(distances))**3/19==F(729,19)
    for p,q in combinations(nodes,2):
        h = squared(p,q)
        if 108<=h<=504:
            i,j = index[p],index[q]
            weight = F(108,h)**3/19
            edges.append((i,j,weight,charge.get((p,q),0)))
            row_sum[i] += weight
            row_sum[j] += weight
            if i<19<=j:
                crossing += 1
                bridge_norm += weight**2
    assert len(nodes)==76 and crossing==114 and max(row_sum)<F(98,1000)
    assert bridge_norm**2/(48*19)>F("2.49e-8")
    print("Crossing edges",crossing,"convexity gap lower bound",float(bridge_norm**2/(48*19)))
    return nodes, edges


def matrices(nodes,edges,w):
    mat = acb_mat(len(nodes),len(nodes))
    for i,j,weight,exponent in edges:
        m = ball(weight)*(w**exponent)
        mat[i,j],mat[j,i] = m,m.conjugate()
    return mat


def trace_values(mat):
    eye = acb_mat([[int(i==j) for j in range(mat.ncols())] for i in range(mat.nrows())])
    value = mat*(eye-mat**3).inv()
    return sum(value[i,i].real for i in range(19)), sum(value[i,i].real for i in range(mat.nrows()))


def calibrated_scale(mat, target):
    # Finite expansion through degree 40 and a geometric bound on the remaining trace.
    cube = mat**3
    power = mat*cube
    coefficients = []
    for j in range(4,41,3):
        coefficients.append(sum(power[i,i].real for i in range(19)))
        if j < 40:
            power = power*cube
    r = ball(F(98,1000))
    remainder = 19*r**43/(1-r**3)
    lo, hi = F(0), F(1)
    for _ in range(100):
        mid = (lo+hi)/2
        t = ball(mid)
        t3 = t**3
        value = arb(0)
        for c in reversed(coefficients):
            value = c+t3*value
        diff = t**4*value+arb(0,remainder)-target
        assert diff > 0 or diff < 0
        if diff > 0:
            hi = mid
        else:
            lo = mid
    return ball(lo).union(ball(hi))


def check_convexity():
    x, y, t = symbols("x y t")
    f = x/(1-x**3)
    assert cancel(diff(f,x,2)-6*x**2*(2+x**3)/(1-x**3)**3)==0
    integral = ((1-t)*(y+t*(x-y))**2).integrate((t,0,1))
    assert cancel(integral-(y+(x-y)/3)**2/2-(x-y)**2/36)==0


def verify():
    ctx.prec = 180
    check_convexity()
    nodes,edges = geometry()
    num, den = fraction(MASS["response"](3))
    roots = []
    for n in range(11485,11491):
        rho = acb.zeta_zero(n)
        w = (rho-1)/rho
        mat = matrices(nodes,edges,w)
        original = acb_mat([[mat[i,j] for j in range(19)] for i in range(19)])
        target,_ = trace_values(original)
        proj,full = trace_values(mat)
        epsilon = MASS["evaluated"](num,w.real)/MASS["evaluated"](den,w.real)
        assert epsilon.overlaps(target)
        assert proj-target>arb(".00003163") and full>arb(".0003976")
        alpha = calibrated_scale(mat,target)
        roots.append(alpha)
        print(n, "projected mass",proj.str(14),"calibrated common factor",alpha.str(26))
    gaps = [u-v for u,v in zip(roots,roots[1:])]
    assert all(g > arb("1.69e-17") for g in gaps)
    print("positive consecutive scale gaps",[g.str(9) for g in gaps])


if __name__ == "__main__":
    verify()
