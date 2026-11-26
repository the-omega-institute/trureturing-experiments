"""用全物理参数域的必要界和区间可行根判定荷六 Hill 截止的最小整数指数。"""
from pathlib import Path
from runpy import run_path
from flint import arb, arb_mat, ctx


def exclude(levels, c, rows, p0, weak_low):
    # 任一合格读出的第 16、18 层已由 K 与 D 定位。
    floor=arb(".14330")
    cap=arb("4.26148e-34")
    alpha=-sum(a for g,a,z in rows)/(2*arb.pi())
    omega=-sum(g for g,a,z in rows)/(8*arb.pi())
    variation=sum(abs(v) for v in c)
    ek=variation*sum(abs(a)*z*(2-z.log()) for g,a,z in rows if not z.is_zero())/(2*arb.pi())
    ew=variation*sum(abs(g)*z for g,a,z in rows)/(4*arb.pi())
    error2=ew/omega+arb("1e-122")/omega+variation*cap**2
    g,b,v=c[16],-c[17],c[18]
    assert min(g,b,v)>0
    e=error2/floor**2
    rmin=((b-g-e)/v).sqrt()
    rmax=((b+e)/v).sqrt()
    # A=x18,B=x19,R=B/A,m=t/(1+t)；以下界对全部 4627≤d≤5020 同时成立。
    mmin=1-1/rmin
    mmax=(1-(1-cap/floor)/rmax)/(1-(levels[17]/levels[18])**4627)
    u17=(levels[16]/levels[17])**4627
    small=cap/floor+u17/mmin
    refined=((b-g*small**2-e)/v).sqrt()
    a_min=(1-ek-alpha*variation*cap)/(alpha*(b-v*refined))
    assert rmin>arb("2.67") and mmin>arb(".62") and mmax<arb(".653")
    assert refined>arb("2.86") and a_min>arb(".221")
    u16max=(levels[15]/levels[17])**4627
    coefficient=(a_min-cap)/(mmax+u16max)
    h_min=p0-(levels[0]/levels[17])**4627/mmin
    position=h_min+coefficient*(levels[15]/levels[17])**5020
    assert coefficient>arb(".338") and position>p0
    # 前 15 权的后缀均为正；按分部求和，弱差余量不超过把它们全放在 p0 的值。
    assert all(sum(c[i:15])>0 for i in range(15))
    ref=arb("-.074256577194663365545")
    tailfloor=(a_min-cap)*(levels[16]/levels[17])**5020
    err=sum(abs(w) for w in c[16:])*(arb("2.400e-21")+arb("7.468e-36")/tailfloor)
    upper=ref+sum(c[:15])*(weak_low(arb(1))-ref)+c[15]*(weak_low(position/p0)-ref)+err
    print("第18层下限、第16层系数下限、第16层/p0 下限与全域弱差上界",
          *[v.str(12) for v in (a_min,coefficient,position/p0,upper)],flush=True)
    assert upper<arb("-3e-6")


def witness(root,levels,c,rows,p0,weak_low):
    heat=run_path(str(root/"ch39-chiral-minimal-neutral-gravity.py"))["expint"]
    moments=run_path(str(root/"ch39-signed-cutoff-gravity.py"))["moments"]
    weak_high=run_path(str(root/"ch39-full-neutral-free-cutoff.py"))["weak_column"]
    mz,mp=arb("91.1876"),arb("1.22089e19")
    p1=(arb("1.176397e14")/mp)**2
    relative=[(v/levels[17])**5021 for v in levels]
    relative[17]=arb(1)
    center=[arb(".636058916"),arb("1.868159790")]
    radius=arb("1e-8")
    box=[v+arb(0,radius.upper()) for v in center]
    def nodes(h,par):
        k,t=par
        for u in relative:
            y=u/(t+u)
            yield h+k*y,y,-k*y/(t+u)
    def jac(h,par):
        out=[[arb(0)]*2 for _ in range(2)]
        for w,(x,dk,dt) in zip(c,nodes(h,par),strict=True):
            factor=[(-z/x).exp() for g,a,z in rows]
            k=sum(a*e for (g,a,z),e in zip(rows,factor,strict=True))/(2*arb.pi())
            v=-x*sum(g*e for (g,a,z),e in zip(rows,factor,strict=True))/(4*arb.pi())
            for i,f in enumerate((k,v)):
                out[i][0]+=w*f*dk
                out[i][1]+=w*f*dt
        return arb_mat(out)
    # 给定两个有理的最低截止底值；不通过任何数值搜索来定义证书。
    qlo,qhi=arb("91.19165"),arb("91.19173")
    hlo,hhi=(qlo/mp)**2,(qhi/mp)**2
    low=hlo.union(hhi)
    assert box[0]>0 and box[0]+low<1 and box[1]>0
    read=[moments(rows,heat,x) for x,_,_ in nodes(low,center)]
    res=arb_mat([[sum(v*w[i] for v,w in zip(c,read,strict=True))-target]
                 for i,target in enumerate((1,arb("1e-122")))])
    pre=jac(low.mid(),center).inv()
    res=pre*res
    defect=arb_mat([[1,0],[0,1]])-pre*jac(low,box)
    row_norm=[sum(abs(t).upper() for t in row) for row in defect.tolist()]
    assert all(abs(z)+v*radius<radius for z,v in zip(res.entries(),row_norm,strict=True))
    def weak(h):
        values=[]
        for x,_,_ in nodes(h,box):
            values.append(weak_low(x/p0) if x<arb("1e-24") else weak_high(x,p0,p1))
        return sum(w*d for w,d in zip(c,values,strict=True))
    dl,dr=weak(hlo),weak(hhi)
    assert dl>0>dr
    qmin=mp*next(nodes(low,box))[0].sqrt()
    s=box[1]*(levels[17]/levels[18])**5021
    print("5021 异号端点、最低能量及 h,κ,s",dl.str(9),dr.str(9),qmin.str(11),
          low.str(12),box[0].str(12),s.str(12),flush=True)
    assert mz<qmin and arb("91.1916")<qmin<arb("91.1918")


def verify(root):
    run_path(str(root/"ch39-charge-six-cutoff-slope.py"))["verify"](root)
    gravity=run_path(str(root/"ch39-three-index-free-atoms.py"))
    _,_,(rows,)=gravity["broken_fields"](root,[2_306_888_100_085])
    native=run_path(str(root/"ch39-native-rh-ray-pvm.py"))
    _,levels,direct,_=next(item for item in native["shell_weights"](root) if item[0]==6)
    ctx.prec=256
    c=[v.real for v in direct[0]]
    p0=(arb("91.1876")/arb("1.22089e19"))**2
    weak=native["weak_series"]()
    exclude(levels,c,rows,p0,weak)
    witness(root,levels,c,rows,p0,weak)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
