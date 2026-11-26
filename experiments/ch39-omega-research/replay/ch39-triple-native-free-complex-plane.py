"""三个实际编号共享完整压缩矩阵的自由复四站三平面区间检验。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, ctx


def adj(m):
    return m.conjugate().transpose()


def col(v):
    return acb_mat([[x] for x in v])


def normalized(v):
    length = (adj(v)*v)[0,0].real
    if not length > 0:
        return None
    return v*(1/length.sqrt())


def measure(y,mat,mass):
    means = [(adj(y)*m*y)[0,0].real.intersection(ms[0].union(ms[-1]))
             for m,ms in zip(mat,mass,strict=True)]
    residuals = [m*y-y*s for m,s in zip(mat,means,strict=True)]
    best = max(range(2),key=lambda q: float((adj(residuals[q])*residuals[q])[0,0].real.lower()))
    z = normalized(residuals[best])
    if z is None:
        return None
    x = col([(y[(i+1)%3,0]*z[(i+2)%3,0]-y[(i+2)%3,0]*z[(i+1)%3,0]).conjugate() for i in range(3)])
    data = []
    for m,ms,s in zip(mat,mass,means,strict=True):
        p = (adj(x)*m*x)[0,0].real.intersection(ms[0].union(ms[-1]))
        w = (adj(z)*m*y)[0,0].real
        k = max(arb(0),abs((adj(x)*m*z)[0,0]).lower())
        if not(p>0 and s>0 and 1-p>0 and 1-s>0 and w.is_finite()):
            return None
        data.append((p,s,w,k*k,ms[0]*ms[1]*ms[2]))
    return data


def line_test(y,mat,mass,ratios):
    data = measure(y,mat,mass)
    if data is None:
        return None
    # 同时翻转 z 的符号等价于交换两个字符与两条正斜率。
    u = max([arb(0)]+[max((w/s).lower(),(-w/(1-s)).lower()) for p,s,w,k2,d in data])
    v = max([arb(0)]+[max((-w/s).lower(),(w/(1-s)).lower()) for p,s,w,k2,d in data])
    low,high = arb(0),1-u*v
    if high < 0:
        return arb(1)
    for (p,s,w,k2,d),r in zip(data,ratios,strict=True):
        g = w/s
        prod = max(arb(0),(u-g).lower())*max(arb(0),(v+g).lower())
        if prod*s*s*p > d:
            return arb(1)
        # b² 的外区间。k2 已经向下取界；其余涉及 p,s 的运算保留完整球。
        lower1 = 4*r*k2/(p*(1-r))**2
        lower2 = k2/(p*(1-p))
        upper1 = d/(s*r*p*p)
        upper2 = 1-w*w/(s*(1-s))
        if not all(t.is_finite() for t in (lower1,lower2,upper1,upper2)):
            return None
        low = max(low,lower1.lower(),lower2.lower())
        high = min(high,upper1.upper(),upper2.upper())
    if low > high or high <= 0:
        return arb(1)
    for (p,s,w,k2,d),r in zip(data,ratios,strict=True):
        sum_bound = low*p*p+2*k2+k2*k2/(high*p*p)
        if r*s*max(4*k2,sum_bound.lower()) > d:
            return arb(1)
    return None


def verify(root):
    scan = run_path(str(root/"ch39-four-site-physical-amplitude-cap.py"))
    ckm = run_path(str(root/"ch39-shell-clock-ckm-blocks.py"))
    ctx.prec = 256
    v = acb_mat(ckm["standard_matrix"](*[arb(x).sin() for x in (".22704",".04275",".00377")],arb("1.139")))
    mass = [list(map(arb,row)) for row in ((".00000615",".00311",".8616"),(".0000135",".000268",".01401"))]
    mat = [adj(v)*scan["diag"](mass[0])*v, scan["diag"](mass[1])]
    normalized_mat = [m*(1/ms[-1]) for m,ms in zip(mat,mass,strict=True)]
    a,b = normalized_mat
    assert not (a*b-b*a).det().contains(0)
    radial = run_path(str(root/"ch39-triple-native-bloch-radius.py"))
    ratios = [arb(value) for value in radial["TRIPLE_RATIO_FLOORS"]]
    for ref,other in [(a,b),(b,a)]:
        pending = [(arb(-1),arb(1),[0,1,2],0)]
        cells,depth = 0,0
        while pending:
            left,right,branches,level = pending.pop()
            keep = scan["pencil_cell"](ref,other,left,right,branches,mat,mass,ratios,line_test=line_test)
            cells += 1
            depth = max(depth,level)
            if keep:
                assert level<45, (left,right,keep)
                middle = (left+right)/2
                pending.extend(((left,middle,keep,level+1),(middle,right,keep,level+1)))
        print("自由复几何一张铅笔图的排除盒数、深度",cells,depth,flush=True)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
