"""保持原谱双引力矩并拟合顶夸克曲率时，三支逆规范耦合的最小极差。"""
from pathlib import Path
from runpy import run_path
import sys
from flint import arb, arb_mat


def verify(root):
    joint = run_path(str(root / "ch39-all-field-signed-crossing.py"))["joint_system"]
    system = joint(root, arb("1e16"), order=16, last=47)
    sites, rows, matrix, target, initial, cols, heat = (
        system[key] for key in ("sites", "rows", "matrix", "target", "initial", "cols", "heat")
    )
    assert not matrix.det().contains(0)
    exact_weights = matrix.solve(target).entries()
    rho = (arb("246.22")/arb("1.22089e19"))**2/2
    y = arb(".8616")
    norm = 1/(8*arb.pi()**2)
    curvature = arb(".258")
    masses = [z for _, _, z in rows[49:65]]
    top = [3*norm*y**4*heat(1,y*y*rho/x) for x in sites]

    def dot(v, w):
        return sum(a*b for a, b in zip(v, w, strict=True))

    def affine(t, direction):
        return [c+t*d for c, d in zip(exact_weights, direction, strict=True)]

    def direction(a, b):
        return matrix.solve(arb_mat([[0], [0], [0], [a], [b]])).entries()

    radius = arb("21.669155")
    unit = arb(-1).union(arb(1))
    bands = [[-norm*heat(1,(z+rho*unit)/x) for x in sites] for z in masses]
    bounds = []
    for a, b in ((1,0), (1,1), (0,1), (-1,0), (-1,-1), (0,-1)):
        weights = affine(radius, direction(a,b))
        lower = dot(weights, top) + sum(min(arb(0), dot(weights, row).lower()) for row in bands)
        assert lower > curvature
        bounds.append(lower.lower())
    print("半径 21.669155 的六顶点曲率下界", *[v.str(13) for v in bounds])

    east = direction(1,0)
    def witness(portals):
        deriv = [-3*norm*y*y*x*heat(2,y*y*rho/x) for x in sites]
        coeff = list(top)
        for lam, z in zip(portals, masses, strict=True):
            for k, x in enumerate(sites):
                deriv[k] += norm*lam*x*heat(2,(z+lam*rho)/x)
                coeff[k] -= norm*lam*lam*heat(1,(z+lam*rho)/x)
        den = dot(east, coeff)
        assert not den.contains(0)
        t = (curvature-dot(exact_weights, coeff))/den
        weights = affine(t,east)
        return t, dot(weights,deriv), weights

    t, _, _ = witness([arb(1)]*12+[arb(0)]*4)
    assert radius < t < arb("21.669158")
    print("只指定曲率的精确见证规范极差", t.str(15))

    def stationary(lam):
        return witness([arb(1)]*7+[-arb(1)]*5+[lam]+[arb(0)]*3)
    lo, hi = arb(".2270635"), arb(".2270639")
    assert stationary(lo)[1] > 0 > stationary(hi)[1]
    t, _, weights = stationary(lo.union(hi))
    assert arb("21.66994") < t < arb("21.66996")
    scores = [a+dot(weights, [c[i+3] for c in cols]) for i, a in enumerate(initial)]
    assert min(a.lower() for a in scores) > 10
    print("同时精确达到驻点和曲率的规范极差、五权、三支逆耦合",
          t.str(15), *[c.str(13) for c in weights], *[a.str(13) for a in scores])


if __name__ == "__main__":
    verify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent)
