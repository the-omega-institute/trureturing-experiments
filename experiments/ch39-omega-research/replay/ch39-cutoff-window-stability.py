"""三枚相邻零点的电弱背景共用截止在六个区间上的必需权重。"""
from functools import cache
from pathlib import Path
from runpy import run_path
from flint import arb,arb_mat


def prepare(root):
    old=run_path(str(root/"ch39-three-index-free-atoms.py"))
    column,fields=old["spectra"](root)
    q=list(map(arb,("1.70333178802223e-9","2.54051739237733e-7","8.96947966979927e-6",
                    ".0303687851337465",".54440151704291","1")))
    v=[column(x) for x in q]
    tangent=[q[1]*z for z in old["column_slope"](fields,q[1])]
    A=arb_mat([*v,tangent])
    F=[v.mid() for v in old["solve"](A,arb_mat([[-1],[0],[-1],[0],[-1],[0],[0]]))]
    L=[]
    for i in range(6):
        L.append([v.mid() for v in old["solve"](A,arb_mat([[int(i==j)] for j in range(7)]))])
    return column,fields,q,F,L,old


def rows_for(fields,d):
    _,ak,aw,bk,bw,ck,cw=d
    rows=[(z,a*ak/(2*arb.pi()),-g*aw/(4*arb.pi())) for g,a,z in fields[0][:4]]
    for field,k,w in zip(fields,(ak-bk-ck,bk,ck),(aw-bw-cw,bw,cw),strict=True):
        rows.extend((z,a*k/(2*arb.pi()),-g*w/(4*arb.pi())) for g,a,z in field[4:])
    return rows


def verify(root):
    column,fields,q,F,L,old=prepare(root)
    api=run_path(str(root/"ch39-common-cutoff-negative-mass.py"))
    tau=arb("1e-6")
    eta=arb("1e-10")
    signs=(-1,1,-1,1,-1,1)
    duals=[[a+tau*s*b for a,b in zip(F,l)] for s,l in zip(signs,L)]
    rows=[rows_for(fields,v) for v in duals]
    big=[api["long_series"](r,22) for r in rows]
    cut=arb("1e-5")
    # 向内取 dyadic 端点，证出的质量也落在声明的闭区间内。
    windows=[((v*arb(".98")).upper(),min(arb(1),(v*arb("1.02")).lower())) for v in q]
    knots=sorted([arb(0),arb("1e-28"),cut,*[v for box in windows for v in box],*q[:-1]],key=float)
    limit=arb("1308843.19305")
    @cache
    def raw(key):
        x=arb(key)
        assert x.is_exact()
        v=column(x)
        return [sum(a*b for a,b in zip(d,v)) for d in duals]
    def evaluate(x):
        # 将端点取成精确 dyadic，使缓存复用同一个函数值。
        return raw(x.mid().man_exp())
    def permitted(v,j,index,margin):
        lo=-1-margin-(tau if j==index and signs[j]<0 else 0)
        hi=margin+(tau if j==index and signs[j]>0 else 0)
        return lo<v<hi
    leaves=0
    for a,b in zip(knots,knots[1:]):
        a,b=a.mid(),b.mid()
        index=next((i for i,(l,r) in enumerate(windows) if l.mid()<=a and b<=r.mid()),None)
        todo=[(a,b,list(range(6)))]
        while todo:
            lo,hi,current=todo.pop()
            lval,hval=evaluate(lo),evaluate(hi)
            interval=lo.union(hi)
            retry=[]
            for j in current:
                assert permitted(lval[j],j,index,eta) and permitted(hval[j],j,index,eta),(j,float(lo),float(hi),lval[j],hval[j])
                der=big[j](interval) if lo>=cut.mid() else api["small_derivative"](interval,rows[j])
                if der>0 or der<0:
                    continue
                enc=lval[j]+der*arb(0).union(hi-lo)
                if not permitted(enc,j,index,eta):
                    retry.append(j)
            if retry:
                assert hi-lo>arb("1e-40")
                m=((lo+hi)/2).mid()
                todo.extend(((lo,m,retry),(m,hi,retry)))
            else:
                leaves+=1
    eps=arb("1e-122")
    goals=[d[0]+d[1]+d[2]*eps for d in duals]
    minmass=[(v-limit-eta*(1+2*limit))/tau for v in goals]
    print("七矩见证截止处允许2%相对误差；全区间认证叶段数",leaves)
    print("负质量不超过",limit.str(17),"时，六个窗口按原子符号所需质量的严格下界",
          *[z.lower().str(16) for z in minmass])
    floors=(1057871,1256561,244466,50019,5691,1449)
    assert all(z>v for z,v in zip(minmass,floors))
    missing=[(val-eta)/(1+2*eta) for val in goals]
    assert min(v.lower() for v in missing)>arb("1308843.19449")
    print("缺任一符号窗口的负量共同下界",min(v.lower() for v in missing).str(20))
    verify_six(column,fields,q,old,limit)


def verify_six(column,fields,q,old,limit):
    # 单列可动的原有六点构造。先缩小切向辅助系数的反号区间，再估计负质量。
    solve=old["solve"]
    cols=[column(v) for v in q]
    tangent=[q[1]*v for v in old["column_slope"](fields,q[1])]
    target=arb_mat([[1],[1],[arb("1e-122")],[0],[0],[0],[0]])
    def weights(v):
        return solve(arb_mat([cols[0],v,tangent,*cols[2:]]).transpose(),target)
    left,right=(q[1]*(1-arb("1e-13"))).mid(),(q[1]*(1+arb("1e-13"))).mid()
    assert weights(column(left))[2]>0>weights(column(right))[2]
    for _ in range(16):
        middle=((left+right)/2).mid()
        value=weights(column(middle))[2]
        if value>0:
            left=middle
        else:
            assert value<0
            right=middle
    middle=((left+right)/2).mid()
    radius=max(middle-left,right-middle).upper()
    box=left.union(right)
    slope=old["column_slope"](fields,box)
    values=weights([v+arb(0,radius)*m for v,m in zip(column(middle),slope)])
    assert all(s*values[j]>0 for j,s in zip((0,1,3,4,5,6),(-1,1,-1,1,-1,1)))
    cost=-values[0]-values[3]-values[5]
    assert cost<limit
    print("六原子构造的严格负量上界",cost.upper().str(22))


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
