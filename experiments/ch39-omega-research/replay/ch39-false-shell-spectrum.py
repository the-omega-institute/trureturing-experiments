"""任意四振幅下，有限错误整壳的谱读数统一下界及两个最近见证的最优拟合。"""
from bisect import bisect_left
from itertools import combinations
from math import factorial, prod
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, arb_mat, ctx


N = 2_306_888_100_085


def spectral_distance(x,y,freq):
    table=[]
    for t in x:
        row=[]
        for w in freq:
            s,c=(t*w).sin_cos()
            row.extend((c,s))
        table.append(row)
    a=arb_mat(table)
    b=arb_mat([[z] for z in y])
    least=(a.transpose()*a).solve(a.transpose()*b)
    residual=b-a*least
    euclidean=sum(z*z for z in residual.entries()).sqrt()
    # 五点上的带符号四行子式是该子表的左零向量。选出一个同时可行的对偶与原问题见证。
    minors={r:arb_mat([table[i] for i in r]).det() for r in combinations(range(len(x)),4)}
    best=None
    for rows in combinations(range(len(x)),5):
        v=[(-1)**i*minors[rows[:i]+rows[i+1:]] for i in range(5)]
        target=sum(z*y[j] for z,j in zip(v,rows))
        lower=abs(target)/sum(abs(z) for z in v)
        if best is None or float(lower)>float(best[0]):
            best=lower,rows,v,target
    lower,rows,v,target=best
    assert all(abs(z)>0 for z in v) and abs(target)>0
    signs=[1 if z*target>0 else -1 for z in v]
    fit=arb_mat([table[i]+[sign] for i,sign in zip(rows,signs)]).solve(arb_mat([[y[i]] for i in rows]))
    assert lower.overlaps(fit[4,0])
    assert all(abs(y[i]-sum(z*fit[k,0] for k,z in enumerate(row)))<lower for i,row in enumerate(table) if i not in rows)
    return lower,euclidean,rows


def distances_to_targets(angles,omegas):
    # 两端都可选作参考。逐个目标作两次有序查找，避免枚举三个弧长的组合。
    result=[]
    for base in angles:
        errors=[]
        for omega in omegas:
            near=[]
            for target in (base-omega,base+omega):
                j=bisect_left(angles,target)
                if j<len(angles):
                    assert angles[j]>target
                    near.append(abs(angles[j]-target))
                if j:
                    assert angles[j-1]<target
                    near.append(abs(angles[j-1]-target))
            errors.append(min(z.lower() for z in near))
        result.append(max(errors))
    return min(result)


def dictionary(root,omegas):
    source=run_path(str(root/"ch39-nineteen-lattice-shells.py"))
    ctx.prec=128
    lower=[arb(1),arb(1)]
    count=0
    dist=source["square_distance"]
    for u in range(1,218,2):
        if u%17==0:
            continue
        pts=source["shell"](3*u*u)
        if len(pts)>217:
            continue
        count+=1
        angles={d:(1-arb(d)/(6*u*u)).acos() for d in {dist(p,q) for p,q in combinations(pts,2)}|{0}}
        lower[0]=min(lower[0],distances_to_targets([angles[d] for d in sorted(angles)],omegas))
        for p in sorted({tuple(sorted(p)) for p in pts}):
            ds=sorted({dist(p,q) for q in pts})
            lower[1]=min(lower[1],distances_to_targets([angles[d] for d in ds],omegas))
    assert count==87
    print("参考弧可取两个方向时的全字典最大角频误差下界",*[v.str(16) for v in lower],flush=True)
    return lower


def verify(root):
    ctx.prec=288
    labels=[(510,1),(1173,-2)]
    theta=[(1-arb(d)/1734).acos() for d in (244,724,804)]
    omega=[z-theta[0] for z in theta[1:]]
    gap=omega[1]-omega[0]
    delta=dictionary(root,omega)
    ctx.prec=288
    h=[z.imag for z in acb.zeta_zeros(N,8)]
    t=[arb(0),*[z-h[0] for z in h[1:]]]
    t.extend(acb.zeta_zero(N+k).imag-h[0] for k in (-279,-232,-45,99))
    c=[]
    for r,i in labels:
        phase=acb(r,i*arb(867).sqrt())
        c.append(-36*(phase/abs(phase)).imag)
    y=[sum(coef*((x+h[0])*w).sin() for coef,w in zip(c,omega)) for x in t]
    assert arb(".27215586483699")<y[0]<arb(".27215586483700")
    x=t[:8]
    lambda2=2*sum(w*w for w in omega)+4*(arb.pi()/2)**2
    time2=sum((v-sum(x)/8)**2 for v in x)
    determinant=1-lambda2*time2/(2*(8**2-1))
    assert determinant>0
    weights=sum(1/prod(abs(a-b) for k,b in enumerate(x) if k!=j) for j,a in enumerate(x))
    print("短窗归一化行列式下界与七阶重心权绝对和",determinant.str(16),weights.str(16))
    for d,floor in zip(delta,("5.76e-18","5.68e-17")):
        assert 0<d<gap/6
        separation=2*omega[0]**3*(omega[0]+omega[1])*gap*(2*gap*d/3)
        lower=min(abs(a).lower()/2 for a in c)*determinant*separation/(factorial(7)*weights)
        assert lower>arb(floor)
        print("有限字典八、十二项逐点谱误差共同严格下界",lower.str(16))
    for u,ds in ((211,(41428,117556,130104)),(105,(196,7756,9576))):
        tip=[(1-arb(d)/(6*u*u)).acos() for d in ds]
        freq=[v-tip[0] for v in tip[1:]]
        print("频差见证壳、平方弦长、两频漂移",u,ds,*[(a-b).str(16) for a,b in zip(freq,omega)])
        for n in (8,12):
            worst,euclidean,rows=spectral_distance(t[:n],y[:n],freq)
            expected={(211,8):"2.031602305897354e-11",(211,12):"1.047112597001384e-5",
                      (105,8):"2.973576009760529e-10",(105,12):"2.675187305014034e-4"}
            assert abs(worst/arb(expected[u,n])-1)<arb("1e-14")
            print("观测数、任意四振幅的最优最大误差、最优欧氏误差、无穷范数活动行",n,
                  worst.str(16),euclidean.str(16),rows,flush=True)


if __name__=="__main__":
    verify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent)
