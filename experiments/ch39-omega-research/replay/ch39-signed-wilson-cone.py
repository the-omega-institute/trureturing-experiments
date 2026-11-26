"""检验共用的有符号 Wilson 面权在六荷正电导下的三迹可行性及负质量代价。"""
from itertools import product
from json import loads
from pathlib import Path
from runpy import run_path
import sys

from flint import arb, arb_mat


def geometry(root):
    old = run_path(str(root/"ch39-wilson-conductance-law.py"))
    _,shell,_,edges,base,faces,models = old["data"](root)
    labels = run_path(str(root/"ch39-mincut-real-triangle-tomography.py"))["enumerate_faces"]()
    sites = shell["geo"]["sites"]
    rotation = [sites.index((y,z,x)) for x,y,z in sites]
    def rotate_edge(e):
        return tuple(sorted(rotation[i] for i in e))
    def orbit(k,turn):
        return tuple(sorted({k,turn(k),turn(turn(k))}))
    position={e:j for j,e in enumerate(edges)}
    face_position={tuple(sorted(s)):j for j,(s,z) in enumerate(faces)}
    def turn_face(i):
        return face_position[tuple(sorted(rotate_edge(e) for e in faces[i][0]))]
    edge_groups=sorted({orbit(i,lambda k:position[rotate_edge(edges[k])]) for i in range(len(edges))})
    groups=sorted({orbit(i,turn_face) for i in range(len(faces))})
    assert len(edge_groups)==45 and len(groups)==179
    edge_class={e:i for i,g in enumerate(edge_groups) for e in g}
    n=len(faces)
    H=arb_mat(3,n)
    terms=[]
    for j,(side,z) in enumerate(faces):
        first,second=side[:2]
        clock,(a,b)=labels[first[0],first[1],second[1]]
        if clock==(0,0) and b==0:
            # 整数相位标签在无钟的实轴上，Wilson 环量精确等于符号。
            z=arb(1 if a>0 else -1)
        values=[1-(z**q).real for q in range(1,7)]
        numbers=[position[e] for e in side]
        weight=2*sum(base[e] for e in side)
        for k,q in enumerate((2,4,6)):
            H[k,j]=weight*values[q-1]
        terms.append((numbers,values))
    Hg=arb_mat([[sum(H[i,j] for j in g) for g in groups] for i in range(3)])
    three=2*sum(base.values())
    fixed=[sum(models[q][3])-three for q in (2,4,6)]
    ranges=[(sum(models[q][2][0])-three,sum(models[q][2][-1])-three) for q in (2,4,6)]
    for q,(lo,hi) in zip((2,4,6),ranges):
        options=models[q][2]
        assert all(lo<sum(v)-three<hi for v in options[1:-1])
    return H,Hg,terms,edge_class,edge_groups,groups,fixed,ranges


def check_dual(record,H,terms,edge_class,edge_groups,targets,cost):
    y=[arb(v)/10**10 for v in record["y"]]
    g=[arb(0)]*270
    for i,numerator in record["mu"]:
        assert numerator>=0
        g[i]=arb(numerator)/10**10
    # 对 529 个独立面分别验不等式，乘子等分给同边轨道的各成员。
    red=[]
    for f,(edges,powers) in enumerate(terms):
        val=sum(H[k,f]*y[k] for k in range(3))
        for edge in edges:
            k=edge_class[edge]
            val+=sum(g[45*q+k]*v/len(edge_groups[k]) for q,v in enumerate(powers))
        red.append(val)
    nonflat=[i for i in range(H.ncols()) if H[0,i]>0]
    flat=set(range(H.ncols()))-set(nonflat)
    assert all(H[0,i]==0 and red[i]==0 for i in flat)
    # 向内移动并正比例缩放，使近似对偶严格满足 -c_f <= y.h_f + mu.B_f <= 0。
    shift=max(arb(0),max((red[i]/H[0,i]).upper() for i in nonflat))+arb("1e-13")
    shifted=[v-shift*H[0,i] for i,v in enumerate(red)]
    weights=[arb(1) if cost=="coefficients" else H[0,i] for i in range(len(red))]
    factor=max(arb(1),max((-shifted[i]/weights[i]).upper() for i in nonflat))+arb("1e-12")
    assert all(-weights[i]<shifted[i]/factor<0 for i in nonflat)
    y[0]-=shift
    return (min(sum(v*t for v,t in zip(y,target)).lower() for target in targets)-sum(g))/factor


def check_witness(record,H,Hg,terms,groups,target):
    a=[arb(0)]*len(groups)
    for i,num in record:
        a[i]=arb(num)/10**7
    # 三个精确定义的修正消除有理近似所留下的三迹误差。
    pivots=(43,60,61)
    sub=arb_mat([[Hg[i,j] for j in pivots] for i in range(3)])
    residual=arb_mat([[v] for v in target])-Hg*arb_mat([[v] for v in a])
    delta=sub.solve(residual).entries()
    for i,v in zip(pivots,delta):
        a[i]+=v
    full=[arb(0)]*H.ncols()
    for group,v in zip(groups,a):
        for j in group:
            full[j]=v
    assert all(v.contains(0) for v in (H*arb_mat([[v] for v in full])-arb_mat([[v] for v in target])).entries())
    ratios=[[arb(1) for _ in range(135)] for _ in range(6)]
    signed=[(-v).max(arb(0)) for v in full]
    for v,(edges,powers) in zip(full,terms):
        for q,power in enumerate(powers):
            for e in edges:
                ratios[q][e]+=v*power
    bottom=min(v.lower() for row in ratios for v in row)
    assert bottom>arb(".0001"),bottom
    mass=sum(v*H[0,j] for j,v in enumerate(signed))
    return sum(signed).upper(),mass.upper(),bottom


def traces_only(H,Hg,terms,edge_class,edge_groups,groups,fixed,corners):
    # 只要求三迹时，每种代价的近最优解均可取两个正面轨道、一个负面轨道。
    for cost,negative in (("coefficients",61),("mass",2)):
        indices=(43,60,negative)
        matrix=arb_mat([[Hg[i,j] for j in indices] for i in range(3)])
        price=arb(len(groups[negative])) if cost=="coefficients" else Hg[0,negative]
        dual=matrix.transpose().solve(arb_mat([[0],[0],[-price]])).entries()
        record={"mu":[],"y":[int((v*10**10).mid().floor().unique_fmpz()) for v in dual]}
        for name,targets in (("fixed",[fixed]),("all",corners)):
            lower=check_dual(record,H,terms,edge_class,edge_groups,targets,cost)
            witness=[]
            for target in targets:
                a,b,c=matrix.solve(arb_mat([[v] for v in target])).entries()
                assert min(a,b)>0 and c<0
                witness.append((-price*c).upper())
            upper=min(witness)
            assert 0 <= (upper-lower).upper() < arb(".02")
            print("仅三迹的负代价、分支、可证下上界",cost,name,lower.str(15),upper.str(15),flush=True)


def verify(root):
    H,Hg,terms,edge_class,edge_groups,groups,fixed,ranges=geometry(root)
    data=loads(Path(__file__).with_suffix(".json").read_text())
    corners=list(product(*ranges))
    # 两个代价的盒下界都于 (低,高,低) 达到其所用对偶线性函数的最小值。
    best=corners[2]
    traces_only(H,Hg,terms,edge_class,edge_groups,groups,fixed,corners)
    results={}
    for cost,col in (("coefficients",0),("mass",1)):
        for kind,values,goal in (("fixed",[fixed],fixed),("box",corners,best)):
            name=cost+"_"+kind
            lower=check_dual(data["edge_dual"][name],H,terms,edge_class,edge_groups,values,cost)
            answer=check_witness(data["witness"][name],H,Hg,terms,groups,goal)
            results[name]=(lower,answer[col])
            print("代价、目标、严格下界、正电导见证上界、最小电导比",cost,kind,
                  lower.str(16),answer[col].str(16),answer[2].str(12),flush=True)
    points=[check_witness(data["witness"]["corner"+str(i)],H,Hg,terms,groups,t) for i,t in enumerate(corners)]
    print("八个目标角点见证最小电导比、最大负面系数和、最大负荷二增迹",
          min(p[2] for p in points).str(14),max(p[0] for p in points).str(14),max(p[1] for p in points).str(14))
    assert min(p[2] for p in points)>arb(".09")
    for threshold,pair in zip(((494310279,494334289),(434171113,434193478),(10813004,10813478),(9510465,9510939)),results.values()):
        assert threshold[0]<pair[0] and pair[1]<threshold[1]


if __name__=="__main__":
    root=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent
    del sys.argv[1:]
    verify(root)
