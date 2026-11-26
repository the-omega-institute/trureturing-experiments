"""最小割三窗内首次允许严格局部极小的正系数三角绕行阶数。"""
from math import cos, sin
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, arb_mat, ctx


def jet(order, row, x, y):
    values = [c*acb(0,order*t).exp() for c,t in zip(row,(x,y,x-y),strict=True)]
    a,b,d = values
    gx,gy = order*(a.imag+d.imag), order*(b.imag-d.imag)
    scale = order**2
    return gx,gy,scale*(a.real+d.real),-scale*d.real,scale*(b.real+d.real)


def prove_instability(model, rows):
    lower, upper = arb("18.94088212"), arb("18.94414547")
    curvature = model["sign"]["density"]["eigen_band"]()
    assert lower < curvature[0] < curvature[1] < upper
    low, high = [model["cell_count"](2,level) for level in (lower,upper)]
    six = model["cell_count"](6,(15-arb(185).sqrt())/2+arb("3e-8"))
    half = arb.pi()/16
    pending = [(1,0,0)]
    excluded = certified = 0
    while pending:
        n,i,j = pending.pop()
        if low(n,i,j)[0]>=16 or high(n,i,j)[1]<16 or six(n,i,j)[1]<1:
            excluded += 1
            continue
        cx,cy = [half*(2*k+1-n)/n for k in (i,j)]
        x,y = cx+arb(0,half/n),cy+arb(0,half/n)
        jets = [jet(q,row,x,y) for q,row in enumerate(rows[:20],1)]
        # 找同一对 u,v，使每一阶的 v·g-uᵀHu > 5|u|²。
        # 每个 u 都取为明确的二进有理向量，三角函数浮点值只影响选择。
        vx,vy = 3*cx/2,3*cy/2
        directions = [(arb(float(cx)*cos(t)-float(cy)*sin(t)),
                       arb(float(cx)*sin(t)+float(cy)*cos(t))) for t in (0,.15,-.15,.3,-.3,.45,-.45)]
        if any(u*u+v*v > 0 and all(
                vx*gx+vy*gy-u*u*hxx-2*u*v*hxy-v*v*hyy > 5*(u*u+v*v)
                for gx,gy,hxx,hxy,hyy in jets) for u,v in directions):
            certified += 1
            continue
        assert n < 2**18, (n,i,j)
        pending.extend((2*n,2*i+a,2*j+b) for a in (0,1) for b in (0,1))
    print("二十阶以下驻点必有曲率小于 -5 的方向，认证格与谱排除格",certified,excluded)


def mincut_observables(model, shell, n):
    ctx.prec = 112
    height = acb.zeta_zero(n).imag
    ctx.prec = 224
    physical = model["sign"]["physical"]
    inputs,initial = physical["model_inputs"]()
    chiral = physical["joint"]["chiral"]
    rows = chiral["source"]["spectral_supplier"](inputs)(height)
    masses = {q:chiral["ordered"]([v for v,k in rows if k==q]) for q in (0,*inputs["charges"])}
    delta = arb(23)/25*arb(2)**-20/(height**2+shell["soft"]["point"]["study"]["offset"])
    cut = 9+(1-delta)**2
    assert all(masses[q][2] < cut < masses[q][3] for q in inputs["charges"])
    joint = physical["joint"]
    _,_,alpha = joint["gauge"]["crossing_enclosure"](inputs,initial,masses,False,724,
                                                   "29.04",".1184",recenter=True)
    vacuum = joint["vacuum"]
    heavy = 32*arb("2.435e18")/724**2*masses[2][15]
    squared = (arb("246.22")**2/(heavy*arb("1.22089e19")))**2
    balanced = vacuum["paired_spectrum"](vacuum["minimal"]["fields"](masses),min_gap="1e-18")
    fields = balanced+vacuum["low_sector"](4,squared)
    assert vacuum["newton"](fields, arb(".03")) < 1 < vacuum["newton"](fields,arb(".033"))
    coefficient = (9*arb(3).log()-12*arb(2).log())/2
    curve = squared**2*coefficient/(arb.pi()*arb("1e-122"))*(1-arb(0).union(squared/(arb(".03")*coefficient)))
    assert arb(".117") < alpha < arb(".119") and arb(".998") < curve < arb("1.002")
    print("第 226358 枚零点的强耦合与曲率比",alpha.str(14),curve.str(14))
    return height,delta


def minimum(model, shell, library, coefficients):
    height,delta = mincut_observables(model,shell,226358)
    phases = shell["geo"]["phases"](height)
    x,y = (phases[k].arg() for k in (724,804))
    rows = library["softened"](coefficients,delta)
    orders = (1,3,21)
    jets = [jet(q,rows[q-1],x,y) for q in orders]
    matrix = arb_mat([[gx,gy,1] for gx,gy,*_ in jets]).transpose()
    weights = matrix.solve(arb_mat([[0],[0],[1]])).entries()
    assert all(w>arb(".03") for w in weights)
    a,b,c = [sum(w*item[j] for w,item in zip(weights,jets)) for j in (2,3,4)]
    gap = ((a-c)**2+4*b*b).sqrt()
    low,high = (a+c-gap)/2,(a+c+gap)/2
    assert low > 13 and high > 50
    # 比较同一固定高度、同一已解出的权，在整个相位环面内与零相位的作用差。
    difference = sum(w*sum(z.real-(z*acb(0,q*t).exp()).real
                     for z,t in zip(rows[q-1],(x,y,x-y),strict=True))
                     for q,w in zip(orders,weights))
    print("三个正权、Hessian 本征值、相对零相位的作用差",
          *[v.str(14) for v in (*weights,low,high,difference)])
    assert difference > 0


def verify(root, library):
    model = run_path(str(root/"ch39-mincut-three-window-angle.py"))
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    ctx.prec = 160
    assert shell["soft"]["point"]["study"]["offset"] > 0 and acb.zeta_zero(1).imag > 14
    cap = arb(23)/25*arb(2)**-20/14**2
    coefficients = library["triangle_coefficients"](shell,21)
    prove_instability(model,library["softened"](coefficients,arb(0).union(cap)))
    minimum(model,shell,library,coefficients)


if __name__ == "__main__":
    root = Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent
    library = Path(sys.argv[2]) if len(sys.argv)>2 else root/"ch39-positive-triangle-charge-action.py"
    verify(root,run_path(str(library)))
