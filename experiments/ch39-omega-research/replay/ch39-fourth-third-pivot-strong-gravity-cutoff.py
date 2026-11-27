"""第四壳第三枢轴同时校准两个引力矩所需的六个大权层截止条件。"""
from pathlib import Path
from runpy import run_path

from flint import arb,ctx


def verify(root, *, number=2_306_888_100_089, certified_region=None):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    _,ordered=lib("adjacent-ordered-charge-six-cutoffs")["verify"](root,refine_node=2_306_888_100_086)[2]
    _,_,pivots,_=next(row for row in lib("native-rh-ray-pvm")["shell_weights"](root,number) if row[0]==6)
    c=[v.real for v in pivots[2]]
    _,heat,(field,)=lib("three-index-free-atoms")["broken_fields"](root,[number])
    ctx.prec=192
    threshold=arb("1e-6")
    strong=[i for i,v in enumerate(c) if abs(v)>4]
    assert [i+1 for i in strong]==[1,10,16,17,18,19]
    other=[i for i in range(19) if i not in strong]
    assert c[3]>0 and all(c[i]<0 for i in other if i!=3)
    b=-c[2]
    bmax=max(-c[i] for i in other if i not in (2,3))
    a=c[3]
    total=sum(v.sqrt() for v in ordered).upper()
    pi=arb.pi()
    slope_k=sum(alpha for _,alpha,_ in field)/(2*pi)
    slope_w=-sum(g for g,_,_ in field)/(8*pi)
    assert slope_k<0<slope_w
    low_tv=sum(abs(c[i]) for i in strong)
    other_tv=sum(abs(c[i]) for i in other)
    # 强权层小于阈值时，直接估计其整个热矩。其余各层保留无质量主项。
    k_err=low_tv*threshold*sum(abs(alpha) for _,alpha,_ in field)/(2*pi)
    w_err=low_tv*threshold**2*sum(abs(g) for g,_,_ in field)/(8*pi)
    xmax=ordered[-1]
    k_err+=other_tv*xmax*sum(abs(alpha)*(1-heat(2,z/xmax)) for _,alpha,z in field)/(2*pi)
    w_err+=other_tv*xmax*xmax*sum(abs(g)*(arb(1)/2-heat(3,z/xmax)) for g,_,z in field)/(4*pi)
    margin=(1-k_err)/abs(slope_k)
    allowance=(arb("1e-122")+w_err)/slope_w
    print("一次矩余量、二次矩容差、下一大负权、总迹",*[v.str(12) for v in (margin,allowance,bmax,total)],flush=True)
    # y、t 是第四层与第三层的振幅。剩余小负权的两个矩记为 A,B。
    # A^3 <= B(Σ |c|sqrt(x))²，且后一个和至多 bmax(T-y-t)。
    stack=[(arb(0),xmax.sqrt().upper(),arb(0),xmax.sqrt().upper())]
    cells=0
    while stack:
        yl,yh,tl,th=stack.pop();cells+=1
        if certified_region is not None and certified_region(yl, tl, total):
            continue
        if yl+tl>total:
            continue
        A=(a*yl**2+margin-b*th**2).lower()
        B=(a*yh**4+allowance-b*tl**4).upper()
        if B<0 or (A>0 and A**3>B*(bmax*(total-yl-tl))**2):
            continue
        assert max(yh-yl,th-tl)>arb("1e-6"),(yl,yh,tl,th)
        if yh-yl>=th-tl:
            mid=(yl+yh)/2
            stack.extend(((yl,mid,tl,th),(mid,yh,tl,th)))
        else:
            mid=(tl+th)/2
            stack.extend(((yl,yh,tl,mid),(yl,yh,mid,th)))
    print("两个例外振幅的闭矩形检查数",cells,flush=True)
    return cells


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
