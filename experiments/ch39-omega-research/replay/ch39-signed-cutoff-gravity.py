"""归一化有符号固有时截止同时匹配牛顿与真空项的负系数总量。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, arb_mat, ctx


NUMBER=2_306_888_100_085


def moments(rows, heat, x):
    if x.is_zero():
        return arb(0),arb(0)
    k=x*sum(a*heat(2,z/x) for g,a,z in rows)/(2*arb.pi())
    w=-x*x*sum(g*heat(3,z/x) for g,a,z in rows)/(4*arb.pi())
    return k,w


def derivatives(rows, x, ell, second=False):
    # 输入区间可在第一小段接触零。各有质量场在零端单侧指数平坦。
    total=arb(0)
    for g,a,z in rows:
        if z.is_zero():
            ex=arb(1)
        elif x.contains(0):
            assert not second
            ex=arb(0).union((-z/x.upper()).exp())
        else:
            ex=(-z/x).exp()
        factor=2*a-ell*g*x
        if second:
            factor=factor*z/(x*x)-ell*g
        total+=ex*factor
    return total/(4*arb.pi())


def certify_sign(rows, ell, lo, hi, sign, second=False):
    pending=[(arb(lo),arb(hi))]
    count=0
    while pending:
        l,r=pending.pop()
        if sign*derivatives(rows,l.union(r),ell,second)>0:
            count+=1
        else:
            assert r-l>arb("1e-21"),(l,r,sign,second)
            m=(l+r)/2
            pending.extend(((l,m),(m,r)))
    return count


def root_box(fn, lower, upper, sign, steps=75):
    l,r=arb(lower),arb(upper)
    assert sign*fn(l)>0>sign*fn(r)
    for _ in range(steps):
        mid=(l+r)/2
        z=sign*fn(mid)
        if z>0:
            l=mid
        elif z<0:
            r=mid
        else:
            break
    return l.union(r)


def weights(rows,heat,points,target):
    cols=[(arb(1),*moments(rows,heat,x)) for x in points]
    mat=arb_mat(cols).transpose()
    assert not mat.det().contains(0)
    return mat.solve(arb_mat([[1],[1],[target]])).entries()


def capped(rows,heat,target, *, return_measure=False):
    a0,a1="1e-6","1.2e-6"
    b0,b1=".49",".51"
    l0,l1=".22150091",".22150096"
    slab=arb(l0).union(arb(l1))
    pieces=sum(certify_sign(rows,slab,lo,hi,sign,sec) for lo,hi,sign,sec in
               (("0",a0,1,False),(a0,a1,-1,True),(a1,b0,-1,False),
                (b0,b1,1,True),(b1,"1",1,False)))
    k1,w1=moments(rows,heat,arb(1))
    def first(ell):
        return root_box(lambda x:derivatives(rows,x,ell),a0,a1,1,62)
    def discrepancy(ell):
        k,w=moments(rows,heat,first(ell))
        return k1-k+ell*(w1-w)
    # 包络定理给差的导数 W(1)-W(a)>0。宽段用积分核 0≤E₃≤1/2。
    assert w1>arb(a1)**2*sum(abs(g) for g,_,_ in rows)/(8*arb.pi())
    ell=root_box(discrepancy,l0,l1,-1,65)
    a=first(ell)
    b=root_box(lambda x:derivatives(rows,x,ell),b0,b1,-1)
    # 直接以极值及两端（一个零点与一个正端点）核对全局支持面。
    ka,wa=moments(rows,heat,a)
    kb,wb=moments(rows,heat,b)
    top,bottom=ka+ell*wa,kb+ell*wb
    assert top>0>bottom and top.overlaps(k1+ell*w1)
    p,n,r=weights(rows,heat,(a,b,arb(1)),target)
    dual=(1+ell*target-top)/(top-bottom)
    assert p>0>n and r>0 and (-n).overlaps(dual)
    print("有界截止：导数覆盖段数、分离斜率及两个内支撑",pieces,ell.str(17),a.str(17),b.str(17))
    print("支撑处的精确解系数与对偶最优负量",*[z.str(17) for z in (p,n,r,dual)])
    if return_measure:
        return dual, (a, b, arb(1)), (p, n, r)
    return dual


def all_phase(root,lib):
    spectrum=run_path(str(root/"ch39-monotone-chiral-phase-obstruction.py"))
    model,*_=spectrum["torus"]["local"]["inputs"]()
    ctx.prec=192
    bands=spectrum["shell_bands"](model)
    wide={q:[v.upper() for v in group] for q,group in bands.items()}
    rows=lib["fields"](wide)
    A=arb(103)/(12*arb.pi())
    B=arb(310)/(8*arb.pi())
    # λ=A/B 时无质量极限为 -Ax(1-x)，全段的质量误差在 x=1 达界。
    heat=lib["expint"]
    err=sum(abs(a)*(1-heat(2,z)) for _,a,z in rows)/(2*arb.pi())
    err+=(A/B)*sum(abs(g)*(arb(1)/2-heat(3,z)) for g,_,z in rows)/(4*arb.pi())
    lower=(1-err)/(A/4+2*err)
    print("全相位包络误差、普朗克有界时统一负量下界",err.str(13),lower.str(13))
    assert lower>arb("1.46316")


def verify(root):
    lib=run_path(str(root/"ch39-chiral-minimal-neutral-gravity.py"))
    chiral=lib["chiral"]
    model,*_=chiral["source"]["mean"]["inputs"]()
    ctx.prec=224
    raw=chiral["source"]["spectral_supplier"](model)(acb.zeta_zero(NUMBER).imag)
    bands={q:chiral["ordered"]([h for h,charge in raw if charge==q]) for q in (0,*model["charges"])}
    rows=lib["fields"](bands)
    assert len(rows)==99 and len(bands[0])==18 and all(len(bands[q])==19 for q in model["charges"])
    assert sum(g for g,_,_ in rows)==-310 and (sum(a for _,a,_ in rows)+arb(103)/6).contains(0)
    assert all(z>0 for _,_,z in rows[1:])
    heat=lib["expint"]
    target=arb("1e-122")
    opt=capped(rows,heat,target)
    assert arb("1.46412119952155")<opt<arb("1.46412119952165")
    for t in (arb(10)**3,arb(10)**6):
        x0=arb("1e-6")
        assert moments(rows,heat,x0)[0]<1
        a,b,c=weights(rows,heat,(x0,t,2*t),target)
        assert a>0>b and c>0
        print("无上限三点例的 t、三个归一系数",t,*[z.str(16) for z in (a,b,c)])
    all_phase(root,lib)


if __name__=="__main__":
    verify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent)
