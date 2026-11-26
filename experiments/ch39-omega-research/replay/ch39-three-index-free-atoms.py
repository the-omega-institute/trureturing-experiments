"""三个相邻零点的破缺引力谱共用可自由选址的有符号截止。"""
from pathlib import Path
from runpy import run_path
from flint import arb, arb_mat, ctx


def broken_fields(root, numbers):
    shared = run_path(str(root/"ch39-two-zero-cutoff-gravity.py"))
    single,lib,original = shared["gravity_fields"](root,numbers)
    rotation = run_path(str(root/"ch39-rotation-portal-top-energy.py"))["rotation_signs"]
    portal = arb(".42995401872")
    r=(arb("246.22")/arb("1.22089e19"))**2/2
    yt,g1,g2=map(arb,(".8616",".467453",".63811"))
    light=[(-59,arb(2)/3,arb(0)),(-12,arb(1),yt**2*r),(6,arb(-1),g2**2*r/2),
           (3,-arb(1)/2,(g2**2+3*g1**2/5)*r/2)]
    total=[]
    for number,field in zip(numbers,original,strict=True):
        signs,masses=rotation(root,number)
        ctx.prec=320
        rows=light+field[1:49]+field[65:]
        for (g,a,z),m,s in zip(field[49:65],masses,signs,strict=True):
            assert z.overlaps(m) and z+s*portal*r>0
            rows.extend(((g-4,a-arb(2)/3,z),(4,arb(2)/3,z+s*portal*r)))
        total.append(rows)
    return single["moments"], lib["expint"], total


def spectra(root):
    n = 2_306_888_100_085
    moments,heat,total = broken_fields(root,(n,n-1,n+1))
    def column(x):
        (k,w),(km,wm),(kp,wp)=[moments(rows,heat,x) for rows in total]
        return [arb(1),k,w,km-k,wm-w,kp-k,wp-w]
    return column,total


def column_slope(fields,x):
    vals=[]
    for rows in fields:
        k,w=arb(0),arb(0)
        for g,a,z in rows:
            exposure=(-z/x).exp()
            k+=a*exposure
            w-=g*x*exposure
        vals.extend((k/(2*arb.pi()),w/(4*arb.pi())))
    k,w,km,wm,kp,wp=vals
    return [0,k,w,km-k,wm-w,kp-k,wp-w]


def solve(a,b):
    n=a.nrows()
    mid=arb_mat([[a[i,j].mid() for j in range(n)] for i in range(n)])
    inverse=mid.inv()
    eye=arb_mat([[int(i==j) for j in range(n)] for i in range(n)])
    e=eye-inverse*a
    norm=max(sum(abs(e[i,j]).upper() for j in range(n)) for i in range(n))
    assert norm<1,norm
    approx=mid.solve(b)
    x=arb_mat([[approx[i,0].mid()] for i in range(n)])
    residual=inverse*(b-a*x)
    radius=max(abs(v).upper() for v in residual.entries())/(1-norm)
    return [v+arb(0,radius.upper()) for v in x.entries()]


def verify(root):
    column,fields=spectra(root)
    q=list(map(arb,("1.70333178802223e-9","2.54051739237733e-7","8.96947966979927e-6",
                    ".0303687851337465",".54440151704291","1")))
    target=arb_mat([[1],[1],[arb("1e-122")],[0],[0],[0],[0]])
    v=[column(x) for x in q]
    h=q[1]*arb("3e-7")
    left,right=q[1]-h,q[1]+h
    def matrix(xcol,spare):
        return arb_mat([v[0],xcol,spare,*v[2:]])
    seven=matrix(column(left),column(right)).transpose()
    c=solve(seven,target)
    assert all(s*w>0 for s,w in zip((-1,1,1,-1,1,-1,1),c,strict=True))
    cost=-c[0]-c[3]-c[5]
    print("七个有理位置的权重及负量",*[w.str(17) for w in c],cost.str(22),flush=True)

    # 固定一列导数，仅用于检验第六个可调位置的单标量余量。
    tangent=[q[1]*s for s in column_slope(fields,q[1])]
    coeff=lambda xx:solve(matrix(xx,tangent).transpose(),target)
    h=q[1]*arb("1e-13")
    left,right=q[1]-h,q[1]+h
    f0,f1=coeff(column(left))[2],coeff(column(right))[2]
    assert f0>0>f1
    dcol=column_slope(fields,left.union(right))
    full=[vv+arb(0,h.upper())*d for vv,d in zip(v[1],dcol,strict=True)]
    box=matrix(full,tangent).transpose()
    six=solve(box,target)
    six_cost=-six[0]-six[3]-six[5]
    assert all(s*six[j]>0 for j,s in zip((0,1,3,4,5,6),(-1,1,-1,1,-1,1),strict=True))
    print("第六原子构造的辅助系数两端及根处负量",f0.str(14),f1.str(14),six_cost.str(18),flush=True)

    interp=matrix(v[1],tangent)
    dual=[d.mid() for d in solve(interp,arb_mat([[-1],[0],[0],[-1],[0],[-1],[0]]))]
    value=(arb_mat([dual])*target)[0,0]
    print("对偶目标与原点值",value.str(18),dual[0].str(16),flush=True)
    _,ak,aw,bk,bw,ck,cw=dual
    rows=[]
    for g,a,z in fields[0][:4]:
        rows.append((z,a*ak/(2*arb.pi()),-g*aw/(4*arb.pi())))
    for field,k,w in zip(fields,(ak-bk-ck,bk,ck),(aw-bw-cw,bw,cw),strict=True):
        rows.extend((z,a*k/(2*arb.pi()),-g*w/(4*arb.pi())) for g,a,z in field[4:])
    api=run_path(str(root/"ch39-common-cutoff-negative-mass.py"))
    small=lambda x:api["small_derivative"](x,rows)
    large=api["long_series"](rows,22)
    edge=arb("1e-5")
    def F(x):
        return sum(d*v for d,v in zip(dual,column(x),strict=True))
    eps=arb("1e-13")
    def permitted(z):
        return z > -1-eps and z < eps
    def cover(a,b,fa,fb,deriv):
        todo=[(a,b,fa,fb)]
        count=0
        while todo:
            l,r,fl,fr=todo.pop()
            assert permitted(fl) and permitted(fr),(l,r,fl,fr)
            interval=deriv(l.union(r))
            if interval>0 or interval<0 or permitted(fl+(arb(0).union(r-l))*interval):
                count+=1
                continue
            assert r-l>arb("1e-40"),(l,r,fl,interval)
            mid=(l+r)/2
            fm=F(mid)
            todo.extend(((l,mid,fl,fm),(mid,r,fm,fr)))
        return count
    limits=[arb(0),arb("1e-28"),*q[:3],edge,*q[3:]]
    values=[F(x) for x in limits]
    for l,r,fl,fr in zip(limits[:-1],limits[1:],values[:-1],values[1:],strict=True):
        count=cover(l,r,fl,fr,small if r<=edge else large)
        print("覆盖区段",l.str(7),r.str(7),count,flush=True)
    lower=(value-eps)/(1+2*eps)
    assert lower>arb("1308843.19304049") and cost<arb("1308843.19304130")
    assert six_cost<arb("1308843.19306")
    print("连续对偶下界、七原子上界",lower.str(18),cost.str(18),flush=True)

if __name__=="__main__":
    import sys
    root=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent
    del sys.argv[1:]
    verify(root)
