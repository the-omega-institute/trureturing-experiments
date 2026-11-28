"""用六质量的正、逆质量形状同时检验 P,H,H² 单个实铅笔的非负谱函数。"""
from itertools import product
from pathlib import Path
from runpy import run_path
import sys

from flint import acb_mat, arb, arb_mat


def shape(values):
    product = values[0]*values[1]*values[2]
    return sum(x*x for x in values)/product**(arb(2)/3)


def mass_shapes(model):
    result = []
    for invert in (False,True):
        for start in (0,3):
            center = model["ref"][start:start+3]
            width = model["tol"][start:start+3]
            corners = [[x+s*w for x,s,w in zip(center,signs,width,strict=True)]
                       for signs in product((-1,1),repeat=3)]
            if invert:
                corners = [[1/x for x in reversed(values)] for values in corners]
            for values in corners:
                # 各坐标盒上最大的质量平方严格支配另外两项。
                assert values[2]**2>2*(values[0]**2+values[1]**2)
            hull = [shape(corner) for corner in corners]
            result.append((min(x.lower() for x in hull),max(x.upper() for x in hull)))
    return result


def eigenbasis(y,z,r):
    length = (1+y*y+z*z).sqrt()
    side = (1+y*y).sqrt()
    vector = arb_mat([[1],[y],[z]])/length
    e = arb_mat([[-y],[1],[0]])/side
    f = arb_mat([[-z],[-y*z],[1+y*y]])/(length*side)
    h,c = y*(1-r*z), z*(r+z)-y*y
    eigenvalue = r+z-r*y*y
    a = eigenvalue-r*h*y-c*z
    matrix = arb_mat([[a,r*h,c],[r*h,r*c,h],[c,h,0]])
    b = 2*(e.transpose()*matrix*f)[0,0]
    d = (e.transpose()*matrix*e-f.transpose()*matrix*f)[0,0]
    gap = (b**2+d**2).sqrt()
    if not gap>0:
        return None
    sign = 1 if d.mid()>=0 else -1
    first,second = gap+sign*d, sign*b
    norm = (first**2+second**2).sqrt()
    if not norm>0:
        return None
    u,w = first/norm, second/norm
    return acb_mat(arb_mat([[vector[i,0],u*e[i,0]+w*f[i,0],-w*e[i,0]+u*f[i,0]]
                            for i in range(3)]))


def coefficients(operators, basis):
    result = []
    for matrix, floor in operators:
        cols = matrix*basis
        norms = [sum(abs(cols[i,j])**2 for i in range(3)) for j in range(3)]
        result.append([(t.lower().max(floor),t.upper()) for t in norms])
    return result


def exclude_weights(coeff, target):
    lo = [max((coeff[k][i][0]/target[k][1]).lower() for k in (2,3)) for i in range(3)]
    hi = [min((target[k][1]/coeff[k][i][0]).upper() for k in (0,1)) for i in range(3)]
    for _ in range(12):
        if any(x>y for x,y in zip(lo,hi,strict=True)):
            return True
        for i in range(3):
            j,k = (i+1)%3,(i+2)%3
            lo[i] = lo[i].max((1/(hi[j]*hi[k])).lower())
            hi[i] = hi[i].min((1/(lo[j]*lo[k])).upper())
        for inverse in (False,True):
            if inverse:
                lo,hi = [(1/x).lower() for x in hi],[(1/x).upper() for x in lo]
            for offset in (0,1):
                row = 2*int(inverse)+offset
                c = coeff[row]
                lower,upper = target[row]
                for i in range(3):
                    j,k = (i+1)%3,(i+2)%3
                    low = (lower-c[j][1]*hi[j]-c[k][1]*hi[k])/c[i][1]
                    high = (upper-c[j][0]*lo[j]-c[k][0]*lo[k])/c[i][0]
                    lo[i],hi[i] = lo[i].max(low.lower()),hi[i].min(high.upper())
                    if lo[i]>hi[i]:
                        return True
            if inverse:
                lo,hi = [(1/x).lower() for x in hi],[(1/x).upper() for x in lo]
    return False


def verify(replay):
    def load(name):
        return run_path(str(replay/("ch39-"+name+".py")))
    model,up,down = load("two-zero-orbit-condition")["calibration"]()
    t,h = load("shell-pencil-six-mass-window")["shell_frame"](model)
    r = h[0,1].real/h[1,2].real
    assert arb(".0006") < r < arb(".0007")
    calibrated = [v*t/abs(v.det())**(arb(1)/3) for v in (up,down)]
    direct = calibrated+[model["adj"](v.inv()) for v in calibrated]
    operators = [(v,(1/sum(abs(x)**2 for x in v.inv().entries())).lower()) for v in direct]
    tar = mass_shapes(model)
    # 权重的迹平均要求至少一条实线的迹比超过下阈值。
    ratio = arb("75.4357")
    assert tar[0][0]/tar[1][1]>ratio
    u,d = [model["adj"](v)*v for v in calibrated]
    q = arb_mat([[(u-ratio*d)[i,j].real for j in range(3)] for i in range(3)])
    block = arb_mat([[-q[i,j] for j in (1,2)] for i in (1,2)])
    assert block[0,0]>0 and block.det()>0
    v = arb_mat([[q[0,1]],[q[0,2]]])
    center = block.solve(v)
    allowance = q[0,0]+(v.transpose()*center)[0,0]
    extent = [(allowance*block.inv()[i,i]).sqrt() for i in range(2)]
    ranges = [(arb(".24"),arb(".56")),(arb(".686"),arb(".719"))]
    y_bound = ranges[0][0].union(ranges[0][1])
    z_bound = ranges[1][0].union(ranges[1][1])
    assert z_bound*(r+z_bound)-y_bound**2>0 and y_bound*(1-r*z_bound)>0
    for i,(left,right) in enumerate(ranges):
        assert left < center[i,0]-extent[i] and center[i,0]+extent[i] < right

    cells = [(*ranges[0],*ranges[1],0)]
    discarded = 0
    covered = 0
    deepest = 0
    while cells:
        ya,yb,za,zb,depth = cells.pop()
        y,z = ya.union(yb),za.union(zb)
        diff = arb_mat([[y-center[0,0]],[z-center[1,0]]])
        if (diff.transpose()*block*diff)[0,0]>allowance:
            discarded += 1
            continue
        frame = eigenbasis(y,z,r)
        if frame is not None and exclude_weights(coefficients(operators,frame),tar):
            covered += 1
            deepest = max(deepest,depth)
            continue
        assert depth<27,(y,z)
        ym,zm = (ya+yb)/2,(za+zb)/2
        cells.extend((a,b,c,d,depth+1) for a,b in ((ya,ym),(ym,yb)) for c,d in ((za,zm),(zm,zb)))
    print("无高迹比线的终格、质量排除终格及最大深度",discarded,covered,deepest)


if __name__=="__main__":
    verify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent)
