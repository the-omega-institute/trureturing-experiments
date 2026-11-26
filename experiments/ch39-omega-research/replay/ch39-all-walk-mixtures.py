"""对全部零、一、两步正混合和任意边相位校验四维一圈投影。"""
from ctypes import CDLL, c_int
from fractions import Fraction as F
from functools import lru_cache
from itertools import combinations
from pathlib import Path
from runpy import run_path
from subprocess import run
from tempfile import TemporaryDirectory

import numpy as np
from flint import acb_mat, arb, arb_poly, ctx


period = 1080
small_period = 270
traces = (F(19), F(23,2), F(1579,216))


def ball(q):
    return arb(q.numerator)/q.denominator


def sparse_nodes():
    cosine, sine, weyl = (np.empty(period) for _ in range(3))
    for i in range(period):
        s,c = (2*arb.pi()*i/period).sin_cos()
        for arr,v in ((sine,s),(cosine,c),(weyl,4*(arb.pi()*i/period).sin()**2)):
            arr[i] = float(v)
            assert abs(v-arb(arr[i])) < arb("1e-15")
            if arr is weyl and i:
                assert abs(v-arb(arr[i])) < abs(v)*arb("1e-15")
    colors = 4*np.array([(a,b,(-a-b)%small_period) for a in range(small_period)
                         for b in range(a+1, small_period) if b<(-a-b)%small_period])
    color_trace = (cosine[colors]+1j*sine[colors]).sum(axis=1)
    wc = np.prod([weyl[(colors[:,i]-colors[:,j])%period]
                  for i,j in ((0,1),(1,2),(0,2))],axis=0)
    weak = np.arange(4,period//2,4)
    angle_rows, weights = [], []
    for y in range(period//6):
        def phase(k):
            return cosine[k%period]+1j*sine[k%period]
        csum = (2*cosine[weak,None]*((phase(y)*color_trace).real+cosine[3*y])
                + ((phase(-4*y)+phase(2*y))*color_trace.conj()).real+cosine[6*y])
        iw, ic = np.nonzero(csum >= 15-2*1.8-1e-10)
        w = weak[iw]
        angles = []
        for c in colors[ic].T:
            angles.extend((y+w+c, y-w+c, -4*y-c, 2*y-c))
        angles.extend((-3*y+w, -3*y-w, np.full_like(w,6*y)))
        angle_rows.append((np.stack(angles,axis=1)%period).astype(np.uint16))
        weights.append(weyl[2*w]*wc[ic])
    nodes, weights = np.concatenate(angle_rows), np.concatenate(weights)
    assert 0 < len(nodes) < 2000000
    return nodes, weights


def vertices(lo,hi):
    points = set()
    for i,t in enumerate(traces):
        if lo <= t <= hi:
            points.add(tuple(F(i == j) for j in range(3)))
    for i,j in combinations(range(3),2):
        for t in (lo,hi):
            p = (t-traces[j])/(traces[i]-traces[j])
            if 0 <= p <= 1:
                point = [F(0)]*3
                point[i],point[j] = p,1-p
                points.add(tuple(point))
    return points


def verify():
    directory = Path(__file__).parent
    shell = run_path(str(directory/"ch39-three-ray-chiral-boundary.py"))
    bounds = run_path(str(directory/"ch39-all-phases-haar-threshold.py"))
    ctx.prec = 384
    upper = bounds["flat_walk"](shell,arb("1e-8"))
    lower = bounds["flat_walk"](shell,arb(0))
    for k,t in enumerate(traces):
        block = lower**k
        assert (sum(block[i,i].real for i in range(19))-ball(t)).contains(0)
    identity = acb_mat(19,19)
    for i in range(19):
        identity[i,i] = 1
    square = upper*upper
    trig = []
    for i in range(period//2+1):
        s,c = (arb.pi()*i/period).sin_cos()
        trig.append((s*s,s*c) if i not in (0,period//2) else (arb(i!=0),arb(0)))

    @lru_cache(maxsize=40)
    def flat_table(point):
        a,b,c = map(ball,point)
        coefficients = (identity*a+upper*b+square*c).charpoly()
        poly = arb_poly([v.real for v in reversed(list(coefficients))])
        derivative = poly.derivative()
        table = []
        for i,(u,du) in enumerate(trig):
            if i == period//2:
                table.append((arb(0),arb(0)))
                continue
            p = poly(u)
            assert p > 0
            table.append((p**6,-du*derivative(u)/p))
        return table

    def build_table(lo,hi):
        flat = [flat_table(p) for p in vertices(lo,hi)]
        rank = traces[1]**2/traces[2]
        max_p0 = (hi-traces[2])/(traces[0]-traces[2])
        variance = (19-rank)/(19*rank)*max(F(0),lo-19*max_p0)**2
        a = ball((lo**2/19+variance)/hi)
        assert 0 < a <= 1
        table = np.empty((period,4))
        for i,(u,du) in enumerate(trig):
            low = min(data[i][0].lower() for data in flat)
            upper_slope = max(data[i][1].upper() for data in flat)
            if i == period//2:
                high = arb(0) if lo == 19 else (1-a)**(6*ball(lo)/a)
                low_slope = arb(0)
            else:
                high = (1-a*u)**(6*ball(lo)/a)
                low_slope = du*ball(lo)/(1-a*u)
            middle = (upper_slope+low_slope)/2
            radius = abs(upper_slope-low_slope)/2
            for j,value in enumerate((low,high,middle,radius)):
                x = float(value)
                assert value == 0 or abs(arb(x)-value) < abs(value)*arb("1e-15")+arb("1e-95")
                table[i,j] = x
                table[-i%period,j] = -x if j == 2 else x
        assert np.max(np.abs(table[:,2:])) < 6600
        return table

    nodes, weights = sparse_nodes()
    small = small_period
    normalization = arb(72)/(period*small**3)
    # Aliasing applies to the actual positive-type face, before bounding it.
    t,u = arb(".15"),arb(".57")
    m_y = sum(d*((q*t).cosh()-1) for q,d in ((1,6),(4,3),(2,3),(3,2),(6,1)))
    tail = 512*((57*m_y-1080*t).exp()+3*(57*8*(u.cosh()-1)-266*u).exp())
    alias = tail*(684**2+5*152**2)
    missing = 256*(-6*ball(traces[2])*arb("1.8")).exp()*(8*19*345)**2
    assert alias+missing < arb("1e-21")
    a3 = 1/arb(".076")
    b1,b2,b3 = arb(41)/10,-arb(19)/6,-7
    # The normalized plane coefficients are separately fractional-linear
    # in sin²(theta_W) and alpha_em^(-1), with denominator of fixed sign.
    for electromagnetic in (arb("127.90"), arb("127.96")):
        for angle in (arb(".23126"), arb(".23132")):
            a1, a2 = 3*(1-angle)*electromagnetic/5, angle*electromagnetic
            assert a1*b2-a2*b1 < 0
            assert (a2*b3-a3*b2)/(a1*b2-a2*b1)>arb(".536")
            assert (a3*b1-a1*b3)/(a1*b2-a2*b1)>-arb("1.516")
    with TemporaryDirectory() as temp:
        lib = Path(temp)/"grid.so"
        run(["clang++","-O3","-std=c++17","-ffp-contract=off","-shared","-fPIC",
             str(directory/"ch39-walk-mixture-grid.cpp"),"-o",str(lib)],check=True)
        project = CDLL(str(lib)).project_nodes
        array = np.ctypeslib.ndpointer
        project.argtypes = [c_int,array(np.uint16,flags="C_CONTIGUOUS"),
                            *[array(np.float64,flags="C_CONTIGUOUS")]*3]
        project.restype = None
        leaves = []
        def cover(lo,hi):
            sums = np.zeros(3)
            project(len(nodes),nodes,weights,build_table(lo,hi),sums)
            rounded = [arb(float(x))*arb("1 +/- 1e-8") for x in sums]
            pos,neg,weak = [x*normalization for x in rounded]
            residual = (pos-neg-alias-missing-arb("1e-75"))/(weak+alias+missing+arb("1e-75"))
            if residual > 0:
                leaves.append((lo,hi,residual.lower()))
                return
            assert hi-lo > F(1,10000), (lo,hi,residual)
            mid = (lo+hi)/2
            cover(lo,mid)
            cover(mid,hi)
        cover(traces[2],traces[0])
        assert all(a[1] == b[0] for a,b in zip(leaves,leaves[1:]))
        assert leaves[0][0] == traces[2] and leaves[-1][1] == traces[0]
        minimum = min(g for _,_,g in leaves)
        assert minimum > arb("3e-5")
        print("完整步数单纯形的迹区间数、节点数、无量纲匹配缺口下界",
              len(leaves),len(nodes),minimum.str(12))


if __name__ == "__main__":
    verify()
