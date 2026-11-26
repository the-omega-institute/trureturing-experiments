"""任意实四站三维子空间的平方质量必要窗与显式正定夸克扩张。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def adj(a):
    return a.conjugate().transpose()


def diagonal(values):
    return acb_mat([[v if i==j else 0 for j in range(len(values))] for i,v in enumerate(values)])


def norm(v):
    return sum(abs(x)**2 for x in v.entries()).sqrt()


def window_bound(mixing):
    # 实际质量窗允许两部门的比例独立缩放。
    masses = [(".00000629",".00316",".8573"), (".0000137",".000271",".01392")]
    radii, heavy_caps = [], []
    for bounds in masses:
        l,m,h = map(arb,bounds)
        u = (m/h)**2
        r = u + (l/h)**2
        radius = (r/(1+r)).sqrt().asin()
        radii.append(radius)
        p = radius.sin()
        heavy_caps.append(u/(1-u-p*p-p))
    total = sum(radii)
    vtd,vts,vtb = [abs(v) for v in mixing[2]]
    angle = vtb.acos()
    assert total < angle and angle+total < arb.pi()/4

    # 两个重方向若同时靠近平面，使用 a<=1 的统一方差条件。
    pt,pb = heavy_caps
    ps = ((vtd+pt+pb)/vts).upper()
    pd = (1-ps**2-pb**2).sqrt().lower()
    required = ((vtd*pd-pt-pb)/vts).lower()
    mean = (arb(".0000137")/arb(".000265"))**2 + ps**2 + (arb(".01410")/arb(".000265"))**2*pb**2
    assert 0 < mean < 1 and required > mean/(1-mean)
    lo,hi = ((angle-total)/2).tan(), ((angle+total)/2).tan()
    print("两条不同固定射线给出 a 的全域必要外界",lo.lower().str(14),hi.upper().str(14))
    return lo.lower(),hi.upper()


def witness(mixing, exponent, geometry, projectors):
    v = acb_mat(mixing)
    masses = [list(map(arb, seq)) for seq in ((".00000615",".00311",".8616"),(".0000135",".000268",".01401"))]
    normalized = [diagonal([(m/seq[-1])**exponent for m in seq]) for seq in masses]
    up,down = adj(v)*normalized[0]*v, normalized[1]
    vals,vecs = (up+down).eig(algorithm="rump",right=True)
    branch = sorted(range(3),key=lambda j:vals[j].real.mid())
    assert all(vals[i].real<vals[j].real for i,j in zip(branch,branch[1:]))
    y = acb_mat([[vecs[i,branch[-1]]] for i in range(3)])
    y = y*(1/norm(y))
    residual = up*y-y*(adj(y)*up*y)[0,0].real
    assert norm(residual)>0
    z = residual*(acb(0,1)/norm(residual))
    x = acb_mat([[(y[(i+1)%3,0]*z[(i+2)%3,0]-y[(i+2)%3,0]*z[(i+1)%3,0]).conjugate()] for i in range(3)])
    unitary = acb_mat([[w[i,0] for w in (x,y,z)] for i in range(3)])
    assert (adj(unitary)*unitary).contains(diagonal([1,1,1]))
    a = geometry
    b = (1-a*a).sqrt()
    e = acb_mat([[1,0,0],[0,0,b],[0,1,0],[0,0,-a]])
    bounds,projected = [],[]
    for target,mass in zip((up,down),masses,strict=True):
        m = adj(unitary)*target*unitary
        assert m[0,1].contains(0) and m[1,2].real.contains(0)
        p,s,t = [m[i,i].real for i in range(3)]
        k,w = m[0,2],m[1,2].imag
        lower = abs(w)/s
        upper = (lower**2+target.det().real/(p*s*s)).sqrt()
        assert lower < a < upper
        bounds.append((lower,upper))
        q = acb_mat([[p,k/b,0,0],[k.conjugate()/b,(t-a*a*s)/(b*b),0,0],
                     [0,0,s,acb(0,-w/a)],[0,0,acb(0,w/a),s]])
        trivial_det = (p*(t-a*a*s)-abs(k)**2)/(b*b)
        assert p>0 and trivial_det>0 and s-abs(w/a)>0
        assert (adj(e)*q*e-m).contains(acb_mat(3,3))
        physical=m*(mass[-1]**exponent)
        sq = physical*physical if exponent==1 else physical
        spectrum,proj=projectors(sq)
        assert all(lam.overlaps(val**2) for lam,val in zip(spectrum,mass))
        projected.append(proj)
    pu,pd=projected
    assert all((pu[i]*pd[j]).trace().real.overlaps(abs(v[i,j])**2) for i in range(3) for j in range(3))
    j=(pu[0]*pd[0]*pu[1]*pd[1]).trace().imag
    assert j>0 and j.overlaps((v[0,0]*v[1,1]*v[0,1].conjugate()*v[1,0].conjugate()).imag)
    low=max(l.upper() for l,h in bounds)
    high=min(h.lower() for l,h in bounds)
    print("目标幂、固定味基的可行 a 内界、有符号 J",exponent,low.str(14),high.str(14),j.str(12))
    return low,high


def verify(root):
    ckm=run_path(str(root/"ch39-shell-clock-ckm-blocks.py"))
    proj=run_path(str(root/"ch39-three-gram-source-flavor.py"))["projections"]
    win=run_path(str(root/"ch39-cyclic-positive-coherent-quarks.py"))["window"]
    ctx.prec=320
    published=(".22704",".04275",".00377","1.139")
    widths=(".00082",".00042",".00008",".023")
    angles=[win(c,w) for c,w in zip(published,widths)]
    literature=ckm["standard_matrix"](*[arb(c).sin() for c in published[:3]],arb(published[3]))
    literature_box=ckm["standard_matrix"](*[c.sin() for c in angles[:3]],angles[3])
    examples = (
        ("书库",ckm["standard_matrix"](*ckm["centers"]),ckm["target_box"](), ".02100345",
         (".0093860",".0326340"), ((".02100338",".02100352"),(".02077540",".03398565"))),
        ("1 TeV 同表",literature,literature_box,".02145737",
         (".0096549",".0332714"), ((".02145730",".02145744"),(".02122434",".03445741"))),
    )
    for name,mixing,box,aval,outer,verified in examples:
        print(name)
        lower,upper = window_bound(box)
        assert arb(outer[0]) < lower < upper < arb(outer[1])
        for power,inner in zip((2,1),verified,strict=True):
            left,right = witness(mixing,power,arb(aval),proj)
            assert left < arb(inner[0]) < arb(aval) < arb(inner[1]) < right


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
