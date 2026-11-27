"""固定十九层实交叉权的保序截止，逐组证实牛顿、真空项及弱三分量差。"""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, ctx


def high_pair(weights,low,rows,heat,moments,*,with_slack=False,scale_text="1e-11"):
    a,b,c=weights
    pi=arb.pi()
    target=arb("1e-122")
    kscale=-sum(t for g,t,z in rows)/(2*pi)
    ratio=(-b/c).sqrt()
    center=[(1/(kscale*(-b-c*ratio))).mid()]
    center.append((ratio*center[0]).mid())
    def values(points,lowx):
        columns=[moments(rows,heat,x) for x in (lowx,*points)]
        return arb_mat([[sum(w*col[j] for w,col in zip(weights,columns,strict=True))-goal]
                        for j,goal in enumerate((1,target))])
    def derivative(points):
        cols=[]
        for w,x in zip((b,c),points,strict=True):
            k=sum(t*(-z/x).exp() for g,t,z in rows)/(2*pi)
            v=-x*sum(g*(-z/x).exp() for g,t,z in rows)/(4*pi)
            cols.append([w*k,w*v])
        return arb_mat(cols).transpose()
    for _ in range(4):
        step=derivative(center).solve(values(center,low.mid()))
        center=[(x-y.mid()).mid() for x,y in zip(center,step.entries(),strict=True)]
    scale=arb(scale_text)
    box=[x+arb(0,scale.upper()) for x in center]
    assert low<box[0]<box[1]<1
    pre=derivative(center).inv()
    defect=arb_mat([[1,0],[0,1]])-pre*derivative(box)
    row_norm=[sum(abs(e).upper() for e in row) for row in defect.tolist()]
    residual=pre*values(center,low)
    assert all(abs(r)+n*scale<scale for r,n in zip(residual.entries(),row_norm,strict=True))
    if with_slack:
        slack = min((scale-abs(r)-n*scale).lower() for r,n in zip(residual.entries(),row_norm,strict=True))
        pre_norm = max(sum(abs(x).upper() for x in row) for row in pre.tolist())
        return box,(slack/pre_norm).lower()
    return box


def exclude_positive_tail(weights,rows):
    # 归一化轮廓在第一截止之前恰取一，在其后受这两个后缀界控制。
    a,b=arb(".912"),arb(".332")
    tails=[sum(weights[i:]) for i in range(1,len(weights))]
    assert all(-a<v<b for v in tails)
    pi=arb.pi()
    k=-sum(t for g,t,z in rows)/(2*pi)
    w=-sum(g for g,t,z in rows)/(8*pi)
    assert k>0 and w>0
    err_k=sum(abs(t)*z*(2-z.log()) for g,t,z in rows if not z.is_zero())/(2*pi)
    err_w=sum(abs(g)*z for g,t,z in rows)/(4*pi)
    assert all(z<1 for g,t,z in rows)
    lin=arb(".97")
    support=(a+b)/(4*lin)+b*(lin-1)
    total=sum(abs(t) for t in weights)
    bound=k*support+total*(err_k+lin*k*err_w/w)
    print("荷二首枢轴的 K+ηW 上界、η、质量修正",bound.str(12),(lin*k/w).str(12),
          err_k.str(10),err_w.str(10),flush=True)
    assert bound<arb(".86")


def verify(root):
    gravity=run_path(str(root/"ch39-three-index-free-atoms.py"))
    _,fields=gravity["spectra"](root)
    rows=fields[0]
    heat=run_path(str(root/"ch39-chiral-minimal-neutral-gravity.py"))["expint"]
    moments=run_path(str(root/"ch39-signed-cutoff-gravity.py"))["moments"]
    native=run_path(str(root/"ch39-native-rh-ray-pvm.py"))
    weak=run_path(str(root/"ch39-full-neutral-free-cutoff.py"))["weak_column"]
    ctx.prec=256
    mz,mp=arb("91.1876"),arb("1.22089e19")
    p0,p1=(mz/mp)**2,(arb("1.176397e14")/mp)**2
    reference=arb("-.074256577194663365545")
    low_weak=native["weak_series"]()
    direct=run_path(str(root/"ch39-neutral-mz-cutoff.py"))["low_weak_column"]
    rho=(arb("246.22")/mz)**2/2
    w=arb(".63811")**2*rho/2
    z=w+3*arb(".467453")**2*rho/10
    t_top=arb(".8616")**2*rho
    splits={(1,3):9,(-3,3):1,(-4,1):17,(6,1):17}
    verified=set()
    for q,eigs,all_weights,_ in native["shell_weights"](root):
        for pivot,raw in enumerate(all_weights,1):
            co=[v.real for v in raw]
            assert sum(co).contains(1)
            key=q,pivot
            if key==(2,1):
                exclude_positive_tail(co,rows)
            if key not in splits:
                continue
            i=splits[key]
            a,b,c=sum(co[:i]),sum(co[i:18]),co[18]
            assert a>1 and b<0<c
            def approx(t):
                return a*low_weak(t)+(1-a)*reference
            left=arb(1)
            right=left
            if q!=-4:
                assert approx(left)>0
                while approx(right)>0:
                    right*=2
                    assert right<arb("1e8")
                assert approx(right)<0
                for _ in range(24):
                    mid=(left+right)/2
                    val=approx(mid)
                    if val>0:
                        left=mid
                    else:
                        assert val<0
                        right=mid
            if q==-4:
                left,right=arb(".89658254"),arb(".89658255")
                low=p0*left.union(1)
            else:
                low=p0*left.union(right)
            mid,high=high_pair((a,b,c),low,rows,heat,moments)
            fixed=b*weak(mid,p0,p1)+c*weak(high,p0,p1)
            if q!=-4:
                dl,dr=[a*low_weak(t)+fixed for t in (left,right)]
                assert dl>0>dr
            else:
                dl,dr=[a*direct(w,z,t_top,p0,p1,x_over_p0=t)+fixed for t in (left,right)]
                assert dl>0>dr
                # D(x) 负且递减。按分部求和，负后缀深度限定可恢复的最大弱差。
                depth=max((-sum(co[k:])).upper() for k in range(1,19))
                upper=low_weak(arb(1))+depth*(low_weak(arb(1))-weak(arb(1),p0,p1))
                print("荷负四首枢轴全保序弱差上界",upper.str(10),flush=True)
                assert upper<arb("-.006")
            print("荷、枢轴、三段长度",q,pivot,i,18-i,1,"组权",*[v.str(12) for v in (a,b,c)],
                  "最低截止 GeV、另两个平方位置",(mz*left.union(right).sqrt()).str(12),
                  mid.str(13),high.str(13),flush=True)
            verified.add(key)
    assert verified==set(splits)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
