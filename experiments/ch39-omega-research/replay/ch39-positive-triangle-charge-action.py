"""最小割三窗内三角作用的力矩，及允许全部绕行数时首次可驻定的阶数。"""
from collections import defaultdict
from itertools import combinations
from math import atan2, cos, prod, sin
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, arb_mat, ctx


def triangle_coefficients(shell, maximum=19):
    edges = {}
    pole = shell["soft"]["pole"]
    tuned = shell["soft"]["star"]["orbit"]
    for i, j, d, z in shell["geo"]["edges"]:
        spoke = pole in (i, j)
        soft = int(spoke and (i in tuned or j in tuned))
        phase = acb(1) if spoke else z
        a, b = int(not spoke and d == 724), int(not spoke and d == 804)
        edges[i,j] = (phase, a, b, soft)
        edges[j,i] = (phase.conjugate(), -a, -b, soft)
    by_cycle = defaultdict(list)
    for i, j, k in combinations(range(19), 3):
        cycle = ((i,j), (j,k), (k,i))
        if all(e in edges for e in cycle):
            zs, aa, bb, ss = zip(*(edges[e] for e in cycle), strict=True)
            a, b, z = sum(aa), sum(bb), prod(zs)
            if a < 0 or (a == 0 and b < 0):
                a, b, z = -a, -b, z.conjugate()
            by_cycle[a,b,sum(ss)].append(z)
    assert {k:len(v) for k,v in by_cycle.items()} == {
        (0,0,0):412, (0,0,1):21, (0,1,0):39, (0,1,1):3,
        (1,0,0):39, (1,0,1):3, (1,-1,0):12}
    coefficients = []
    for q in range(1, maximum+1):
        row = []
        for a, b in ((1,0), (0,1), (1,-1)):
            soft = sum((z**q for z in by_cycle[a,b,1]), acb(0))
            total = soft+sum(z**q for z in by_cycle[a,b,0])
            row.append((total,soft))
        coefficients.append(row)
    return coefficients


def softened(coefficients, delta):
    return [[total-delta*soft for total,soft in row] for row in coefficients]


def gradient(q, row, x, y):
    def term(c, t):
        return c.real*(q*t).sin()+c.imag*(q*t).cos()
    mixed = term(row[2], x-y)
    return q*(term(row[0],x)+mixed), q*(term(row[1],y)-mixed)


def verify_radial(coefficients):
    half = arb.pi()/16
    radius = arb(".10124")
    minima = []
    for q, (a,b,d) in enumerate(coefficients[:12], 1):
        assert 0 < q*arb(".256255") < arb.pi() and d.real > 0
        drift = q*((abs(a.imag).upper()+abs(d.imag).upper())**2+(abs(b.imag).upper()+abs(d.imag).upper())**2).sqrt()
        force = q*min(a.real.lower(),b.real.lower())*radius*(q*half).sin()/half-drift
        assert force > arb("3.00565")
        minima.append(force)
    print("荷 1、2、3、4、6 的向外径向力下界", *[minima[q-1].lower().str(9) for q in (1,2,3,4,6)])


def verify_cover(model, coefficients):
    # 三窗必要条件。按现有荷二曲率谱带和荷六最小割谱带逐格排除。
    lower, upper = arb("18.94088212"), arb("18.94414547")
    curvature = model["sign"]["density"]["eigen_band"]()
    assert lower < curvature[0] < curvature[1] < upper
    low, high = [model["cell_count"](2, level) for level in (lower, upper)]
    six = model["cell_count"](6, (15-arb(185).sqrt())/2+arb("3e-8"))
    floor = arb("3.00565")
    half = arb.pi()/16
    pending = [(1,0,0)]
    excluded = accepted = 0
    while pending:
        n, i, j = pending.pop()
        if low(n,i,j)[0] >= 16 or high(n,i,j)[1] < 16 or six(n,i,j)[1] < 1:
            excluded += 1
            continue
        cx,cy = [half*(2*k+1-n)/n for k in (i,j)]
        radius = half/n
        x,y = cx+arb(0,radius),cy+arb(0,radius)
        if x**2+y**2 > arb(".10124")**2:
            vectors = [gradient(q,row,x,y) for q,row in enumerate(coefficients[:18],1)]
            angle = atan2(float(cy),float(cx))
            # 浮点中点只选择候选方向；每次接受必须通过整个闭格的区间不等式。
            directions = [(round(4096*cos(angle+j/16)),round(4096*sin(angle+j/16)))
                          for j in range(-24,25)]
            centers = [(float(dx),float(dy)) for dx,dy in vectors]
            a,b = max(directions,key=lambda v:min(v[0]*dx+v[1]*dy for dx,dy in centers))
            if all(a*dx+b*dy > floor*arb(a*a+b*b).sqrt() for dx,dy in vectors):
                accepted += 1
                continue
        assert n < 2**20, (n,i,j)
        pending.extend((2*n,2*i+a,2*j+b) for a in (0,1) for b in (0,1))
    print("全体 1 至 18 次绕行的共同投影下界、合格外包格与排除格",floor.str(8),accepted,excluded)


def stationary(model, shell, coefficients):
    n = 295872666901
    ctx.prec = 112
    height = acb.zeta_zero(n).imag
    ctx.prec = 224
    physical = model["sign"]["physical"]
    inputs, initial = physical["model_inputs"]()
    chiral = physical["joint"]["chiral"]
    masses = chiral["source"]["spectral_supplier"](inputs)(height)
    bands = {q:chiral["ordered"]([v for v,k in masses if k == q]) for q in (0,*inputs["charges"])}
    assert all(bands[q][2] < 9 and bands[q][3] > 10 for q in inputs["charges"])
    _, alpha, curve = physical["layer_observables"](inputs,initial,bands)
    assert arb(".117") < alpha < arb(".119") and arb(".998") < curve < arb("1.002")
    delta = arb(23)/25*arb(2)**-20/(height**2+shell["soft"]["point"]["study"]["offset"])
    phases = shell["geo"]["phases"](height)
    x,y = (phases[k].arg() for k in (724,804))
    rows = softened(coefficients,delta)
    orders = (2,3,19)
    vectors = [gradient(q,rows[q-1],x,y) for q in orders]
    matrix = arb_mat([[dx,dy,1] for dx,dy in vectors]).transpose()
    weights = matrix.solve(arb_mat([[0],[0],[1]])).entries()
    assert all(w>arb(".16") for w in weights)
    hessian = arb_mat(2,2)
    for weight,q in zip(weights,orders,strict=True):
        a,b,d = rows[q-1]
        ca,cb,cd = [(c*acb(0,q*t).exp()).real for c,t in zip((a,b,d),(x,y,x-y),strict=True)]
        hessian += weight*q*q*arb_mat([[ca+cd,-cd],[-cd,cb+cd]])
    assert hessian.det() < -370000
    print("第 295872666901 枚的三窗、相角与驻点正权", alpha.str(9),curve.str(10),x.str(10),y.str(10),
          *[w.str(14) for w in weights])
    print("驻点的 Hessian 行列式",hessian.det().str(15))


def verify(root):
    model = run_path(str(root/"ch39-mincut-three-window-angle.py"))
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    ctx.prec = 160
    assert shell["soft"]["point"]["study"]["offset"] > 0 and acb.zeta_zero(1).imag > 14
    delta_cap = arb(23)/25*arb(2)**-20/14**2
    coefficients = triangle_coefficients(shell)
    uniform = softened(coefficients,arb(0).union(delta_cap))
    verify_radial(uniform)
    verify_cover(model,uniform)
    stationary(model,shell,coefficients)


if __name__ == "__main__":
    verify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent)
