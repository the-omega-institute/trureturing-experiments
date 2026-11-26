"""原十九点 Wilson 面权在六种正边电导下保持奇荷迹的检验。"""
from itertools import product
from pathlib import Path
from runpy import run_path
import sys

from flint import arb, arb_mat


def geometry(root):
    baseline=run_path(str(root/"ch39-signed-wilson-cone.py"))
    H,Hg,terms,edge_class,edges,groups,fixed,ranges=baseline["geometry"](root)
    native=run_path(str(root/"ch39-wilson-conductance-law.py"))["data"](root)
    _,shell,_,ed,base,faces,_=native
    B=arb_mat(270,len(groups))
    for j,g in enumerate(groups):
        for f in g:
            incident,powers=terms[f]
            for e in incident:
                i=edge_class[e]
                if e==edges[i][0]:
                    for q,v in enumerate(powers):
                        B[45*q+i,j]+=v
    total=2*sum(base.values())
    W=[2*len(g)*base[ed[g[0]]] for g in edges]
    return B,W,total,fixed,ranges,edges,groups,terms,H,ed,base,shell


def sub(A, rows, columns):
    return arb_mat([[A[i,j] for j in columns] for i in rows])


def column(values):
    return arb_mat([[v] for v in values])


def certify_coordinates(root, data, B, edges, groups, terms, ed, shell):
    labels=run_path(str(root/"ch39-mincut-real-triangle-tomography.py"))["enumerate_faces"]()
    def key(f):
        xyz=tuple(sorted({v for e in terms[f][0] for v in ed[e]}))
        clock,(a,b)=labels[xyz]
        return clock,a,abs(b) if clock==(0,0) else b
    # 轮换保持每个面的全部实 Wilson 幂，亦保持边的软辐条类别。
    assert all(len({key(f) for f in g})==1 for g in groups)
    pole=shell["soft"]["pole"]
    soft=shell["soft"]["star"]["orbit"]
    assert all(len({pole in ed[e] and any(v in ed[e] for v in soft) for e in g})==1 for g in edges)
    C,R=data["C"],data["R"]
    assert sorted(set(range(179))-set(C))==[111,171]
    assert len(C)==len(set(C))==len(R)==len(set(R))==177
    assert key(groups[111][0])==((0,0),1,0)
    # 每组有相同的相位，且带号入射差都沿第 35 与 40 边轨道的差方向。
    bundles=[
        [(139,1),(141,-1)], [(158,1),(159,-1)],
        [(170,1),(156,-1),(169,-1),(154,1)], [(164,-1),(175,1)],
        [(173,1),(151,-1),(148,1),(171,-1)], [(142,-1),(143,1)],
        [(167,1),(145,1),(166,-1),(149,-1)]]
    assert len({j for row in bundles for j,k in row})==sum(map(len,bundles))
    for combo in bundles:
        assert len({key(groups[j][0]) for j,k in combo})==1
        h=[sum(k*sum(g[0] in terms[f][0] for f in groups[j]) for j,k in combo) for g in edges]
        assert set(i for i,v in enumerate(h) if v)=={35,40} and abs(h[35])==1 and h[40]==-h[35]
    # 其中无第 171 列的六组独立，故第七组也在它们张成的六维空间内。
    inverse=sub(B,R,C).inv()
    return sub(B,range(270),C)*inverse, inverse


def trace_rows(U,W,t0):
    return arb_mat([[sum(w*U[45*q+i,j] for i,w in enumerate(W))/t0
                     for j in range(177)] for q in range(6)])


def moment(weights,ed):
    degrees=[arb(0) for _ in range(19)]
    for (i,j),v in zip(ed,weights):
        degrees[i]+=v
        degrees[j]+=v
    return sum(v*v for v in degrees)+2*sum(v*v for v in weights)


def positive_witness(mode, record, scale, J, K, target, t0, U, inv, C, groups, edges, terms, ed, base):
    x=column([arb(v)/scale for v in record])
    b=column([arb(0)]*(K.nrows()-3)+[t/t0 for t in target])
    # 精确见证 x 定义为有理种子加上以下可逆小系统的解。
    x=pin_linear(x,K,J,b)
    assert all(v.contains(0) for v in (K*x-b).entries())
    # 展开为共享的 529 个原面系数并直接计算全部 810 个电导。
    coeff=(inv*x).entries()
    full=[arb(0)]*529
    for j,a in zip(C,coeff):
        for f in groups[j]:
            full[f]=a
    r=[[arb(1) for _ in ed] for _ in range(6)]
    for a,(incident,powers) in zip(full,terms):
        for e in incident:
            for q,p in enumerate(powers):
                r[q][e]+=a*p
    # 个别面权的区间宽度含相关性，用已证明等价的主元坐标收紧每条实际边的区间。
    change=(U*x).entries()
    for q in range(6):
        for i,orbit in enumerate(edges):
            for e in orbit:
                r[q][e]=r[q][e].intersection(1+change[45*q+i])
    floor=min(v.lower() for row in r for v in row)
    assert floor>arb(".99" if mode=="e1" else ".000009"),floor
    native=[base[e] for e in ed]
    weights=[[w*z for w,z in zip(native,row)] for row in r]
    T=[2*sum(row)/t0 for row in weights]
    assert all((T[q]-1-t/t0).contains(0) for q,t in zip((1,3,5),target))
    if mode=="e1":
        assert all((v-1).contains(0) for v in r[0])
    else:
        assert all((T[q]-1).contains(0) for q in (0,2))
    Q=[moment(weights[q],ed)/moment(native,ed) for q in (0,2,4)]
    return floor,T,Q


def certify_lower(rec,scale,U,K,L,targets):
    y=arb_mat([[arb(v)/scale for v in rec["y"]]])
    mu=arb_mat(1,270)
    for i,z in rec["mu"]:
        assert z>0
        mu[0,i]=arb(z)/scale
    residual=L-y*K-mu*U
    norm=sum(abs(v).upper() for v in residual.entries())
    bound=min((y*column(b))[0,0].lower() for b in targets)
    return 1+bound-sum(mu.entries())-500*norm,norm


def pin_linear(x,K,J,b):
    correction=sub(K,range(K.nrows()),J).solve(b-K*x)
    for j,v in zip(J,correction.entries()):
        x[j,0]+=v
    return x


def certify_four(record,scale,J,K,U,L,points,t0,edges,ed,base):
    # 先在一组选定轻根下，沿五个迹等式的核寻找同时满足两个二次矩的点。
    b=column([0,0,*[v/t0 for v in points[record["corner"]]]])
    x=pin_linear(column([arb(v)/scale for v in record["seed"]]),K,J,b)
    basis=[]
    for j in record["directions"]:
        direction=arb_mat(177,1)
        direction[j,0]=1
        basis.append(pin_linear(direction,K,J,column([0]*5)))
    reach=[(U*v).entries() for v in [x,*basis]]
    original=[base[e] for e in ed]
    q0=moment(original,ed)
    def bilinear(u,v):
        # 极化可沿用相同的二次矩计算。
        return (moment([a+b for a,b in zip(u,v)],ed)-moment(u,ed)-moment(v,ed))/(2*q0)
    constant=[]
    linear=[]
    quad=[]
    for q in (0,2):
        vectors=[]
        for j,a in enumerate(reach):
            weights=[arb(0)]*135
            for i,g in enumerate(edges):
                for e in g:
                    weights[e]=original[e]*(int(j==0)+a[45*q+i])
            vectors.append(weights)
        A,P,Q=vectors
        constant.append(moment(A,ed)/q0-1)
        linear.append([2*bilinear(A,P),2*bilinear(A,Q)])
        quad.append([moment(P,ed)/q0,2*bilinear(P,Q),moment(Q,ed)/q0])
    # 精确的二次系统 c+Dt+H(t,t)=0。连续映射 -D^{-1}(c+H(t,t)) 把闭盒送入内部。
    radius=arb("1e-7")
    interval=arb(0,radius.upper())
    close=-(arb_mat(linear).inv()*(column(constant)+arb_mat(quad)*column([interval*interval]*3)))
    assert all(abs(v)<radius for v in close.entries()),close
    anchor=x+basis[0]*close[0,0]+basis[1]*close[1,0]
    # 在每一个目标角点逐边固定同一组荷一、荷三电导，故插值也保持这两个二次矩。
    O=arb_mat(sub(U,[*range(45),*range(90,135)],range(177)).tolist()+sub(K,(2,3,4),range(177)).tolist())
    traits=sub(O,range(90),range(177))*anchor
    lower=arb(1)
    fifth=[]
    for seed,goal in zip(record["transport"],points):
        rhs=column([*traits.entries(),*[v/t0 for v in goal]])
        z=pin_linear(column([arb(v)/scale for v in seed]),O,record["J"],rhs)
        minimum=min((1+v).lower() for v in (U*z).entries())
        assert minimum>arb(".05428"),minimum
        lower=min(minimum,lower)
        fifth.append(1+(sub(L,[4],range(177))*z)[0,0])
    assert min(v.lower() for v in fifth)>arb("10.51252") and max(v.upper() for v in fifth)<arb("11.38481")
    print("荷一、三的一二阶同时保持，二维不动点盒内像上界",max(abs(v).upper() for v in close.entries()).str(10),
          "全部角点共有的最小边比下界",lower.str(10),"第五荷一阶比范围",
          min(v.lower() for v in fifth).str(12),max(v.upper() for v in fifth).str(12),flush=True)


def verify(root):
    from json import loads
    data=loads(Path(__file__).with_suffix(".json").read_text())
    B,W,t0,fixed,ranges,edges,groups,terms,H,ed,base,shell=geometry(root)
    assert shell["model"]["charges"]==(1,4,2,3,6)
    assert len(ed)==135 and len(terms)==529
    U,inv=certify_coordinates(root,data,B,edges,groups,terms,ed,shell)
    L=trace_rows(U,W,t0)
    K=sub(L,(0,2,1,3,5),range(177))
    E=arb_mat(sub(U,range(45),range(177)).tolist()+sub(K,(2,3,4),range(177)).tolist())
    points=[fixed,*product(*ranges)]
    scale=data["scale"]
    examples={}
    for mode,eq in (("e1",E),("c13",K)):
        rows=[]
        for n,goal in enumerate(points):
            answer=positive_witness(mode,data[mode][n],scale,data["J"][mode],eq,goal,t0,U,inv,data["C"],groups,edges,terms,ed,base)
            floor,T,Q=answer
            print(mode,n,"最小边比",floor.str(10),"奇荷一阶比",*[T[i].str(12) for i in (0,2,4)],
                  "奇荷二阶比",*[z.str(12) for z in Q],flush=True)
            rows.append(answer)
        examples[mode]=rows
    certify_four(data["four"],scale,data["J"]["c13"],K,U,L,points[1:],t0,edges,ed,base)
    # 若 T5/T0≤10，六荷的所有代表边比均小于 10 T0/min W，因而各主元增量的绝对值小于 500。
    assert max(t[1]/t0+1 for t in ranges)<10 and max(t/t0+1 for t in fixed)<10
    assert 10*t0/min(w.lower() for w in W)<500
    q0=moment([base[e] for e in ed],ed)
    for case,goals,row,lower,upper,second in (
        ("fixed",[fixed],0,"8.80851160","8.80854735","74.27724"),
        ("box",points[1:],1,"8.16269502","8.16273069","63.78489")):
        targets=[[0,0,*[v/t0 for v in t]] for t in goals]
        lb,norm=certify_lower(data["dual"][case],scale,U,K,sub(L,[4],range(177)),targets)
        ub=examples["c13"][row][1][4]
        assert arb(lower)<lb<ub<arb(upper)<10
        # Cauchy 分别应用于 19 个度和 135 个电导。
        second_lb=(arb(1)/19+arb(1)/270)*(lb*t0)**2/q0
        assert second_lb>arb(second)
        print(case,"第五荷一阶下界",lb.str(16),"正见证上界",ub.str(16),
              "对偶残差一范数上界",norm.str(9),"第五荷二阶比下界",second_lb.str(16),flush=True)


if __name__=="__main__":
    here=Path(__file__).resolve().parent
    root=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else here
    del sys.argv[1:]
    verify(root)
