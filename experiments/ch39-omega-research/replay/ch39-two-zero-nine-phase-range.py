"""用正逆相对谱迹的隐函数确定六质量三角度全窗的可达 CP 相位区间。"""
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import arb


def mixing_squares(angles, cosine):
    a,b,z = [v.sin() for v in angles]
    ca,cb,cz = [v.cos() for v in angles]
    cross = 2*a*b*z*ca*cb*cosine
    return ((ca**2*cz**2, a**2*cz**2, z**2),
            (a**2*cb**2+ca**2*b**2*z**2+cross, ca**2*cb**2+a**2*b**2*z**2-cross, b**2*cz**2),
            (a**2*b**2+ca**2*cb**2*z**2-cross, ca**2*b**2+a**2*cb**2*z**2+cross, cb**2*cz**2))


def trace_minors(values, cosine, traces, dual):
    """在满足两条迹等式的点上计算 D_{+,i}D_{-,0}-D_{-,i}D_{+,0}。"""
    weights = [[(x**2/(xs[(i+1)%3]*xs[(i+2)%3]))**(arb(2)*s/3) for i,x in enumerate(xs)]
               for xs,s in ((values[:3],1), (values[3:6],-1))]
    up,down = weights
    ij = tuple(product(range(3),repeat=2))
    ratio = {(i,j):up[i]*down[j] for i,j in ij}
    probs = mixing_squares(values[6:],cosine)
    deriv = []
    for axis in range(6):
        entries = [(i,j) for i,j in ij if (i==axis if axis<3 else j==axis-3)]
        # 对数质量导数中，几何均值的贡献可直接写为已知的目标迹。
        deriv.append([2*s*(1 if axis<3 else -1)*(
            sum(probs[i][j]*ratio[i,j]**s for i,j in entries)-target/3)
            for s,target in zip((1,-1),traces,strict=True)])
    for axis in range(3):
        dm = mixing_squares([dual(x,int(i==axis)) for i,x in enumerate(values[6:])],cosine)
        deriv.append([sum(ratio[i,j]**s*dm[i][j].d for i,j in ij) for s in (1,-1)])
    p0,n0 = deriv[0]
    return [p*n0-n*p0 for p,n in deriv[1:]]


def verify():
    replay = Path(__file__).resolve().parent
    library = run_path(str(replay/"ch39-nine-window-cp-cosine.py"))
    trace = library["affine_trace"]
    dual = library["Dual"]
    normalized_traces = run_path(str(replay/"ch39-two-zero-sharp-quark-fit.py"))["normalized_traces"]
    model,up,down = run_path(str(replay/"ch39-two-zero-orbit-condition.py"))["calibration"]()
    background = normalized_traces(up,down)
    centers,widths = model["ref"][:9],model["tol"][:9]
    signs = (-1,1,-1,1,1,1,1,-1)

    def values(coords):
        return [c+w*x for c,w,x in zip(centers,widths,coords,strict=True)]

    def residual(v, targets):
        a,b = trace(v)
        ai,bi = trace(v,inverse=True)
        cosine = (-a+targets[0])/b
        return cosine, ai+bi*cosine-targets[1]

    unit = arb(-1).union(arb(1))
    for x in (-3,3):
        assert -x*residual(values([x]+[unit]*8),background)[1] > 0

    def segment(lo,hi,i,n):
        left = lo+arb(hi-lo)*i/n
        return left.union(lo+arb(hi-lo)*(i+1)/n)

    for k,l,h in product(range(8),range(4),range(4)):
        coords = [unit]*9
        coords[0],coords[3],coords[8] = segment(-3,3,k,8),segment(-1,1,l,4),segment(-1,1,h,4)
        v = values(coords)
        cosine,f = residual(v,background)
        if not f.contains(0):
            continue
        gradient = residual([dual(t,int(i==0)) for i,t in enumerate(v)],background)[1].d
        assert gradient < 0 and trace(v)[1] < 0
        minors = trace_minors(v,cosine,background,dual)
        assert all(-s*m > 0 for s,m in zip(signs,minors,strict=True))

    for i in range(32):
        t = segment(-1,1,i,32)
        for x in (-1,1):
            assert -x*residual(values([x]+[s*t for s in signs]),background)[1] > 0

    at_zero = run_path(str(replay/"ch39-quadratic-gram-source.py"))["at_zero"]
    actual = normalized_traces(at_zero(model,-4,2345666),at_zero(model,2,964890))
    assert all(b.contains(a) for b,a in zip(background,actual,strict=True))
    for label,targets,brackets,phase_bounds in (
        ("原制备盒", background, (("-.3542","-.3532"),(".0584",".0593")),
         (("1.0413","1.0415"),("1.2362","1.2364"))),
        ("实际零点", actual, (("-.353721123208048","-.353721123208047"),
                             (".0588550999407690",".0588550999407700")),
         (("1.04136193991962","1.04136193991964"),("1.23626570405977","1.23626570405979")))):
        endpoints = []
        for pm,(lo,hi),(minimum,maximum) in zip((1,-1),brackets,phase_bounds,strict=True):
            corner = [pm*s for s in signs]
            lo,hi = arb(lo),arb(hi)
            assert residual(values([lo]+corner),targets)[1] > 0
            assert residual(values([hi]+corner),targets)[1] < 0
            cosine = residual(values([lo.union(hi)]+corner),targets)[0]
            assert -1 < cosine < 1
            phase = cosine.acos()
            assert arb(minimum) < phase < arb(maximum)
            endpoints.append(phase.str(20))
        print(label, "可达绝对相位的最小值与最大值", *endpoints)


if __name__ == "__main__":
    verify()
