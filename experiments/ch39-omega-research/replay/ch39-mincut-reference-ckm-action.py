"""最小割三窗和参考 CKM 同时精确时，绕数三角作用的局部极小阶数。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, arb_mat, ctx


def reference_angles(flavor):
    # 六面中的两种钟边走向相反，故平均实部为 A cos(x-y)+B sin(x-y)。
    def coordinate(z):
        return flavor["singlet_at"]({724: z, 804: acb(1)})[1]

    a, b = coordinate(acb(1)), coordinate(acb(0, 1))
    reference = flavor["singlet"](11661)[1]
    offset = acb(a, b).arg()
    angle = (reference / (a*a + b*b).sqrt()).acos()
    assert arb(".14965576381076") < angle < arb(".14965576381078")
    return offset - angle, offset + angle


def spectral_cut(root, level):
    inertia = run_path(str(root / "ch39-neutral-mean-phase-inertia.py"))
    fixed, terms = inertia["schur_coefficients"](2,level)

    def count(x,y,gap):
        u,v,w = [acb(0,2*t).exp() for t in (x,y,-gap)]
        ranges = [c+2*(a*u+b*v+d*w).real for c,a,b,d in terms]
        low = 16-fixed+sum(t<0 for t in ranges)
        high = 16-fixed+sum(not(t>0) for t in ranges)
        return low,high

    return count


def excludes_pair(a,b):
    ax,ay,axx,axy,ayy = a
    bx,by,bxx,bxy,byy = b
    if not (ax*by-ay*bx).contains(0):
        return True
    if ax*bx > 0 or ay*by > 0:
        return True
    for u,v in ((ax,bx),(ay,by)):
        if (u-v).contains(0):
            continue
        weight = -v/(u-v)
        if weight < 0 or weight > 1:
            return True
        xx,xy,yy = [weight*first+(1-weight)*second
                    for first,second in ((axx,bxx),(axy,bxy),(ayy,byy))]
        if xx < 0 or yy < 0 or xx*yy-xy*xy < 0:
            return True
    return False


def branch_cover(root, source, tri, offsets):
    cap = arb(23)/25*arb(2)**-20/14**2
    coeffs = [source["coeff"](tri,q,arb(0).union(cap)) for q in range(1,32)]
    cuts = [arb("18.94088212"),arb("18.94414547")]
    required = run_path(str(root/"ch39-curvature-window-density-cap.py"))["eigen_band"]()
    assert cuts[0] < required[0] < required[1] < cuts[1]
    low_count, high_count = [spectral_cut(root,t) for t in cuts]
    groups = (arb("-.142").union(arb("-.1365")),arb("-.0132").union(arb("-.0082")))
    for gap in offsets:
        sign = 1 if gap>0 else -1
        positive = gap*sign
        pending = [(1,0)]
        leaves = [0,0,0]
        while pending:
            n,i = pending.pop()
            # y 为带符号的新坐标；实际相位为 sign*(x,y)。
            y = arb.pi()/16*(arb(2*i+1,1)/n-1)
            x = y+positive
            if x > arb.pi()/16 or low_count(sign*x,sign*y,gap)[0]>=16 or high_count(sign*x,sign*y,gap)[1]<16:
                leaves[0]+=1
                continue
            if groups[0].contains(y):
                axis = 1
                active = (1,29,30)
            elif groups[1].contains(y):
                axis = 0
                active = (1,28,29)
            else:
                axis = None
            if axis is not None:
                rows, values = [], []
                for q in active:
                    gx,gy,hxx,_,hyy = source["derivatives"]([(q, coeffs[q-1])],
                                                             sign*x.mid(), sign*y.mid())
                    rows.append([gx,gy,arb(1)])
                    values.append([hyy if axis else hxx])
                u,v,_ = arb_mat(rows).solve(arb_mat(values)).entries()
                u,v = u.mid(),v.mid()
                verified=True
                jets=[]
                for q,cs in enumerate(coeffs,1):
                    a,b,d=[c*acb(0,q*t).exp() for c,t in zip(cs,(sign*x,sign*y,gap),strict=True)]
                    gx,gy=q*(a+d).imag,q*(b-d).imag
                    second=q*q*((b if axis else a)+d).real
                    jets.append((gx,gy,q*q*(a+d).real,-q*q*d.real,q*q*(b+d).real))
                    if q<=30 and not second-(u*gx+v*gy) < -30:
                        verified=False
                        break
                if verified and all(excludes_pair(first,jets[-1]) for first in jets[:-1]):
                    leaves[axis+1]+=1
                    continue
            assert n < 2**21, (sign,n,i,x,y)
            pending.extend([(2*n,2*i),(2*n,2*i+1)])
        assert min(leaves)>0
        print("参考角差、谱排除格与两个坐标负曲率证书的格数",gap.str(15),leaves)


def full_windows(physical,x,y):
    model, initial = physical["model_inputs"]()
    ctx.prec = 224
    joint = physical["joint"]
    supply = joint["chiral"]["source"]["spectral_supplier"](model)
    rows = supply(arb(10)**12,clocks={724:acb(0,x).exp(),804:acb(0,y).exp()})
    bands = {q: joint["chiral"]["ordered"]([v for v,charge in rows if q==charge])
             for q in (0,*model["charges"])}
    assert all(bands[q][2]<arb("9.99999999") and bands[q][3]>10 for q in model["charges"])
    _,_,alpha=joint["gauge"]["crossing_enclosure"](
        model, initial, bands, False, 724, "29.04", ".1184", recenter=True)
    vac = joint["vacuum"]
    heavy = 32*arb("2.435e18")/724**2*bands[2][15]
    ratio = (arb("246.22")**2/(heavy*arb("1.22089e19")))**2
    paired = vac["paired_spectrum"](vac["minimal"]["fields"](bands),min_gap="1e-18")
    fields = paired+vac["low_sector"](4,ratio)
    assert vac["newton"](fields,arb(".03"))<1<vac["newton"](fields,arb(".033"))
    factor = (9*arb(3).log()-12*arb(2).log())/2
    curve = ratio**2*factor/(arb.pi()*arb("1e-122"))*(1-arb(0).union(ratio/(arb(".03")*factor)))
    assert arb(".117") < alpha < arb(".119") and arb(".998") < curve < arb("1.002")
    print("取定相位的强耦合和曲率比",alpha.str(14),curve.str(16))


def verify(root, winding):
    source = run_path(str(winding))
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    physical = run_path(str(root/"ch39-first-three-mode-observables.py"))
    flavor = physical["joint"]["flavor"]
    ctx.prec = 192
    tri = source["terms"](shell)
    gaps = reference_angles(flavor)
    ref = flavor["singlet"](11661)[1]
    branch_cover(root,source,tri,gaps)
    y = arb("-.139")
    x = y+gaps[1]
    state,exact = flavor["singlet_at"]({724:acb(0,x).exp(),804:acb(0,y).exp()})
    assert exact.overlaps(ref)
    cp = run_path(str(root/"ch39-cp-ensemble-test.py"))
    matrix = flavor["reference_matrix"]()
    p,m,j = cp["matrix_readings"](matrix,state)
    p0,m0,j0 = cp["matrix_readings"](matrix,flavor["singlet"](11661)[0])
    assert all(u.overlaps(v) for u,v in zip([*p,*m,j],[*p0,*m0,j0],strict=True))
    full_windows(physical,x,y)
    # 选择实际第 20000 枚的刚度参数。两只联络角独立变分。
    h = acb.zeta_zero(20000).imag
    assert h > 10000
    offset = shell["soft"]["point"]["study"]["offset"]
    delta = arb(23)/25*arb(2)**-20/(h*h+offset)
    orders = (2,30,31)
    coeffs = [source["coeff"](tri,q,delta) for q in orders]
    ders = [source["derivatives"]([(q,c)],x,y) for q,c in zip(orders,coeffs,strict=True)]
    basis = arb_mat([[1]*3,[v[0] for v in ders],[v[1] for v in ders]])
    weights = basis.solve(arb_mat([[1],[0],[0]])).entries()
    assert min(w.lower() for w in weights)>arb(".18")
    hs = [sum(w*v[i] for w,v in zip(weights,ders,strict=True)) for i in (2,3,4)]
    xx,xy,yy = hs
    assert xx>7 and (xx-7)*(yy-7)>xy*xy
    comparison = sum(w*source["energy"](q,c,x,y) for q,c,w in zip(orders,coeffs,weights,strict=True))
    assert comparison > arb("1.97854")
    print("取定正权、Hessian的最小特征值及相对零相位能量",*[w.str(17) for w in weights],
          ((xx+yy-((xx-yy)**2+4*xy**2).sqrt())/2).str(16),comparison.str(16))


if __name__=="__main__":
    root = Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent
    winding = Path(sys.argv[2]) if len(sys.argv)>2 else root/"ch39-mincut-global-winding.py"
    verify(root,winding)
