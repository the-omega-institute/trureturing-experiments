"""公共规范谱平滑对上下夸克质量比的有限前缀检验。"""
from pathlib import Path
from runpy import run_path

from itertools import combinations

from flint import acb, arb, ctx


local = run_path(str(Path(__file__).with_name("ch39-local-light-higgs-hierarchy.py")))


def spectra(q, height, center):
    # 对易的轮换在每条简单本征线上为三次单位根。下游只用互异轮换荷与本征值顺序。
    shell = local["local"]["shell"]
    sites = shell["geo"]["sites"]
    rotate = [sites.index((b,c,a)) for a,b,c in sites]
    lap = shell["laplacian"](q,height)
    assert all((lap[rotate[i],rotate[j]]-lap[i,j]).contains(0) for i in range(19) for j in range(19))
    eig, vec = lap.eig(right=True,algorithm="rump")
    assert all(not z.overlaps(w) for j,z in enumerate(eig) for w in eig[:j])
    omega = acb(-1,arb(3).sqrt())/2
    rows=[]
    for j,z in enumerate(eig):
        mass=abs(z.real*(z.real-center))
        assert mass<9 or mass>9
        if mass<9:
            norm=sum(abs(vec[i,j])**2 for i in range(19)).sqrt()
            col=[vec[i,j]/norm for i in range(19)]
            char=sum(v.conjugate()*col[rotate[i]] for i,v in enumerate(col))
            labels=[k for k in range(3) if char.overlaps(omega**k)]
            assert len(labels)==1
            rows.append((z.real,labels[0],col))
    rows.sort(key=lambda r:float(r[0]))
    assert all(rows[j][0]<rows[j+1][0] for j in range(2))
    assert {r[1] for r in rows}=={0,1,2}
    return rows


def energies(left,right,H,conjugate):
    entries={}
    for i,(_,a,u) in enumerate(left):
        for j,(_,b,v) in enumerate(right):
            for _,k,phi in H:
                mode=[z.conjugate() for z in phi] if conjugate else phi
                val=sum(x*y*h for x,y,h in zip(u,v,mode))
                allowed=(a+b+(-k if conjugate else k))%3==0
                if allowed:
                    entries[i,j,k]=abs(val)**2
                    assert entries[i,j,k]>0
                else:
                    assert val.contains(0)
                    entries[i,j,k]=arb(0)
    return entries


def ordered_prefixes(terms,r):
    ranked=sorted(terms,key=lambda row:float(row[1]*r+row[2]))
    values=[t[1]*r+t[2] for t in ranked]
    assert all(a<b for a,b in zip(values,values[1:]))
    U,D=arb(0),arb(0)
    heights=[]
    for sector,_,_,w in ranked:
        if sector=="u":
            U+=w
        else:
            D+=w
        if D.is_zero():
            return None
        heights.append((U/D).sqrt().upper())
    return max(heights)


def relative_time(left, up, down, u, d, k):
    terms=[(tag, lv, rv, w[i,j,k]) for tag,right,w in (("u",up,u),("d",down,d))
           for i,(lv,_,_) in enumerate(left) for j,(rv,_,_) in enumerate(right) if w[i,j,k]>0]
    assert len(terms)==6
    breaks=[]
    for a,b in combinations(terms,2):
        l1,r1=a[1:3]
        l2,r2=b[1:3]
        if l1.overlaps(l2):
            assert l1 is l2
            continue
        hit=(r2-r1)/(l1-l2)
        if hit>0:
            breaks.append(hit)
        else:
            assert hit<0
    breaks.sort(key=float)
    assert all(a<b for a,b in zip(breaks,breaks[1:]))
    closed=[arb(0),*breaks]
    ranges=[]
    unlimited=[]
    for j,low in enumerate(closed):
        r=(low+closed[j+1]).mid()/2 if j+1<len(closed) else low.upper()+1
        cap=ordered_prefixes(terms,r)
        if cap is None:
            unlimited.append(low)
        else:
            ranges.append(cap)
    return max(ranges),unlimited


def compare(number):
    height=acb.zeta_zero(number).imag
    L0=local["local"]["shell"]["laplacian"](0,height)
    center=sum(L0[i,i].real for i in range(19))/18
    left,up,down,H=[spectra(q,height,center) for q in (1,-4,2,3)]
    u,d=energies(left,up,H,False),energies(left,down,H,True)
    order=sorted([(row[0],q,j) for q,rows in (("u",up),("d",down))
                  for j,row in enumerate(rows)],key=lambda r:float(r[0]))
    assert all(p[0]<q[0] for p,q in zip(order,order[1:]))
    # 右谱共同非增倍率的平方是六个初段指示函数的非负线性组合。
    extrema=[]
    for end in range(1,7):
        active=order[:end]
        for k in range(3):
            U=sum(u[i,j,k] for _,sector,j in active if sector=="u" for i in range(3))
            D=sum(d[i,j,k] for _,sector,j in active if sector=="d" for i in range(3))
            assert D>0
            extrema.append((U/D if U else arb(0)).sqrt())
    maximum=max(x.upper() for x in extrema)
    # 对左谱再施任意正倍率，其条件数乘上述界。κ_L 是振幅倍率之比。
    target=arb("61.36519")
    required=target/maximum
    print(number, "右侧任意单调平滑后的最大范数比、左侧必要条件数",
          maximum.str(10), required.str(10))
    caps,critical=[],[]
    for k in range(3):
        cap,hits=relative_time(left,up,down,u,d,k)
        caps.append(cap)
        critical.extend(hits)
    first=min(critical,key=float)
    largest=max(caps)
    print("联合谱单调响应：仍有有限界的各区间中最大上限、首个无界区间的左/右尺度比",
          largest.str(9),first.str(9))
    assert maximum<arb("1.557") and required>39 and largest<target
    return maximum,required


if __name__=="__main__":
    ctx.prec=384
    for n in (22027,28325):
        compare(n)
