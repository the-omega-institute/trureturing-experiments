"""从四个算术三平面到可匹配夸克数据的实、复四站超平面的主角代价。"""
from pathlib import Path
from runpy import run_path
from flint import arb, ctx


def complex_outer(mixing):
    # A<=4 时，各生成线在横向 y 上的相对斜率至多为 4。
    A = arb(4)
    u,c,t,d,s,b = map(arb,(".00000629",".00316",".8573",".0000137",".000271",".01392"))
    mass2 = [(u/t)**2,(c/t)**2,(d/b)**2,(s/b)**2]
    ru,rd = [(mass2[i]+mass2[i+1]).sqrt().atan() for i in (0,2)]
    vtd,vts,vtb = map(abs,mixing[2])
    theta = vtb.acos()
    assert ru+rd < theta and theta+ru+rd < (1/A).atan()
    caps=[]
    for m,r in zip((mass2[1],mass2[3]),(ru,rd),strict=True):
        den=1-m-r.sin()**2-A*r.sin()
        assert den>0
        caps.append((A*m/den).upper())
    pt,pb=caps
    cd,cs,cb = map(abs,mixing[1])
    ps = arb(1)
    for _ in range(8):
        ps = min(ps,((vtd+vtb*pb+pt)/vts).upper())
        pc = min(arb(1),(cd+cs*ps+cb*pb).upper())
        mu = mass2[0]+mass2[1]*pc**2+pt**2
        md = mass2[2]+mass2[3]*ps**2+pb**2
        assert mu<1 and md<1
        pt=min(pt,(A*mu/(1-mu)).upper())
        pb=min(pb,(A*md/(1-md)).upper())
        ms=(arb(".0000137")/arb(".000265"))**2+ps**2+(arb(".01410")/arb(".000265"))**2*pb**2
        if ms<1:
            ps=min(ps,(A*ms/(1-ms)).upper())
    needed=(vtd*(1-ps**2-pb**2).sqrt()-vtb*pb-pt)/vts
    assert needed>arb(".19") and ps<arb(".03")
    return ((theta+ru+rd)/2).tan().upper()


def certify_unbalanced(a0, angle):
    # 数学区间 [lo,hi] 上，(4+a) 增，sqrt((1-a²)/(16+a²)) 减。
    b0=(1-a0*a0).sqrt()
    pending=[(arb(0),arb(1))]
    while pending:
        lo,hi=pending.pop()
        maximum=a0*hi+b0*(4+hi)*((1-lo*lo)/(2*(16+lo*lo))).sqrt()
        if maximum < angle.cos():
            continue
        assert hi-lo>arb("1e-10")
        mid=(lo+hi)/2
        pending.extend(((lo,mid),(mid,hi)))


def verify(root):
    ckm=run_path(str(root/"ch39-shell-clock-ckm-blocks.py"))
    gram=run_path(str(root/"ch39-native-rh-ray-pvm.py"))["gram_cholesky"](root)
    normal=run_path(str(root/"ch39-four-gram-triplet-ckm.py"))["normal"]
    free=run_path(str(root/"ch39-free-four-site-positive-extensions.py"))
    window=run_path(str(root/"ch39-cyclic-positive-coherent-quarks.py"))["window"]
    ctx.prec=320
    table=[".22704",".04275",".00377","1.139"]
    widths=[".00082",".00042",".00008",".023"]
    angles=[window(x,y) for x,y in zip(table,widths,strict=True)]
    tablebox=ckm["standard_matrix"](*[x.sin() for x in angles[:3]],angles[3])
    bookbox=ckm["target_box"]()
    projections=run_path(str(root/"ch39-three-gram-source-flavor.py"))["projections"]
    book=ckm["standard_matrix"](*ckm["centers"])
    literature=ckm["standard_matrix"](*[arb(v).sin() for v in table[:3]],arb(table[3]))
    data=(("书库",bookbox,book,".0326340",".02100352"),
          ("同尺度表 2",tablebox,literature,".0332714",".02145744"))
    for label,box,center,limit,edge in data:
        left,right=free["witness"](center,2,arb(edge),projections)
        assert left<arb(edge)<right
        outer=arb(limit)
        _, checked=free["window_bound"](box)
        cap=complex_outer(box)
        assert checked<outer and cap<outer
        print(label,"平衡或受控复法向的 a 上界",cap.str(12))
        for omit,bad in enumerate((".52619",".51578",".50435",".49205")):
            n=[v.real for v in normal(gram,omit)]
            a0=(n[0]**2+sum(n[1:])**2/3).sqrt()
            lower=a0.asin()-outer.asin()
            upper=a0.asin()-arb(edge).asin()
            escape=arb(bad)
            certify_unbalanced(a0,escape)
            complex_lower=min(lower.lower(),escape)
            print("省略源、平凡分量、实下界、中心可行上界、复下界",omit,
                  a0.str(12),lower.lower().str(12),upper.upper().str(12),complex_lower.str(12))
            assert (lower>arb(".51338") and complex_lower>arb(".49204") and upper<arb(".595"))


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
