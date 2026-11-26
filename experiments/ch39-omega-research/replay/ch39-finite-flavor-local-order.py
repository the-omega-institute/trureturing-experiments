"""确定六面坐标截断码排除三十圈局部真空所需的最少位数。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, arb_mat, ctx


def certify_cell(source, previous, coefficients, x, y, gap, axis):
    active = (1,29,30) if axis else (1,28,29)
    rows, rhs = [], []
    for q in active:
        gx,gy,hxx,_,hyy = source["derivatives"]([(q,coefficients[q-1])],x.mid(),y.mid())
        rows.append([gx,gy,1])
        rhs.append([hyy if axis else hxx])
    try:
        u,v,_ = arb_mat(rows).solve(arb_mat(rhs)).entries()
    except ZeroDivisionError:
        return False
    u,v = u.mid(),v.mid()
    jets=[]
    for q,cs in enumerate(coefficients,1):
        a,b,d=[c*acb(0,q*t).exp() for c,t in zip(cs,(x,y,gap),strict=True)]
        gx,gy=q*(a+d).imag,q*(b-d).imag
        hxx,hxy,hyy=q*q*(a+d).real,-q*q*d.real,q*q*(b+d).real
        if q<=30 and not (hyy if axis else hxx)-u*gx-v*gy < -1:
            return False
        jets.append((gx,gy,hxx,hxy,hyy))
    return all(previous["excludes_pair"](first,jets[-1]) for first in jets[:-1])


def certify_strips(root, source, previous, tri, phase, amplitude):
    lower,upper = (arb(i)/512 for i in (505,506))
    alpha,beta = (upper/amplitude).acos(),(lower/amplitude).acos()
    cuts = [arb("18.94088212"),arb("18.94414547")]
    required = run_path(str(root/"ch39-curvature-window-density-cap.py"))["eigen_band"]()
    assert cuts[0] < required[0] < required[1] < cuts[1]
    low,high = [previous["spectral_cut"](root,t) for t in cuts]
    delta = arb(0).union(arb(23)/25*arb(2)**-20/14**2)
    coeffs = [source["coeff"](tri,q,delta) for q in range(1,32)]
    for sign in (-1,1):
        # d=sign*(x-y) 是正角距。与 y 的矩形一起用精确二分覆盖取闭包的码格。
        full_d = (sign*phase+alpha).union(sign*phase+beta)
        assert arb(".145329") < full_d < arb(".158262")
        pending = [(full_d.lower(),full_d.upper(),-arb.pi()/16,arb.pi()/16)]
        leaves=[0,0,0]
        while pending:
            dl,dr,yl,yr = pending.pop()
            d,y = dl.union(dr),yl.union(yr)
            x=sign*(y+d)
            y=sign*y
            gap=sign*d
            if abs(x)>arb.pi()/16 or low(x,y,gap)[0]>=16 or high(x,y,gap)[1]<16:
                leaves[0]+=1
                continue
            axis = 1 if abs(y.mid())>abs(x.mid()) else 0
            if max(dr-dl,yr-yl)<arb(".002") and certify_cell(source,previous,coeffs,x,y,gap,axis):
                leaves[axis+1]+=1
                continue
            assert max(dr-dl,yr-yl)>arb("1e-8"), (x,y)
            if dr-dl>yr-yl:
                mid=(dl+dr)/2
                pending.extend(((dl,mid,yl,yr),(mid,dr,yl,yr)))
            else:
                mid=(yl+yr)/2
                pending.extend(((dl,dr,yl,mid),(dl,dr,mid,yr)))
        print("角差符号、曲率和二阶排除叶格", sign,leaves,flush=True)


def verify(root, winding, reference):
    source, previous = run_path(str(winding)),run_path(str(reference))
    physical=run_path(str(root/"ch39-first-three-mode-observables.py"))
    shell=run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    flavor=physical["joint"]["flavor"]
    ctx.prec=192
    tri=source["terms"](shell)
    a,b=[flavor["singlet_at"]({724:z,804:acb(1)})[1] for z in (acb(1),acb(0,1))]
    ref=flavor["singlet"](11661)[1]
    assert (512*ref).floor().unique_fmpz()==505 and (256*ref).floor().unique_fmpz()==252
    certify_strips(root,source,previous,tri,acb(a,b).arg(),abs(acb(a,b)))
    # 八位的反例以一个真实高度固定软辐边，两只联络角自由选取。
    height=acb.zeta_zero(20000).imag
    assert height>10000
    offset=shell["soft"]["point"]["study"]["offset"]
    delta=arb(23)/25*arb(2)**-20/(height**2+offset)
    x,y=arb(".023"),arb("-.143")
    c=flavor["singlet_at"]({724:acb(0,x).exp(),804:acb(0,y).exp()})[1]
    assert (256*c).floor().unique_fmpz()==252 and (512*c).floor().unique_fmpz()==504
    previous["full_windows"](physical,x,y)
    orders=(4,29,30)
    jets=[source["derivatives"]([(q,source["coeff"](tri,q,delta))],x,y) for q in orders]
    matrix=arb_mat([[1]*3,[v[0] for v in jets],[v[1] for v in jets]])
    weights=matrix.solve(arb_mat([[1],[0],[0]])).entries()
    assert min(v.lower() for v in weights)>arb(".23")
    xx,xy,yy=[sum(w*v[k] for w,v in zip(weights,jets,strict=True)) for k in (2,3,4)]
    assert xx>185 and (xx-185)*(yy-185)>xy**2
    print("八位同码见证的坐标、权重、Hessian 最小特征值",c.str(16),
          *[v.str(16) for v in weights],((xx+yy-((xx-yy)**2+4*xy**2).sqrt())/2).str(16))


if __name__=="__main__":
    root=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent
    winding=Path(sys.argv[2]) if len(sys.argv)>2 else root/"ch39-mincut-global-winding.py"
    reference=Path(sys.argv[3]) if len(sys.argv)>3 else root/"ch39-mincut-reference-ckm-action.py"
    verify(root,winding,reference)
