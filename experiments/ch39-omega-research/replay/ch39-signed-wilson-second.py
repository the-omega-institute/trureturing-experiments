"""八组三迹校正后的 Wilson 面权覆盖全部轻根选择的二阶矩锥。"""
from itertools import product
from json import loads
from pathlib import Path
from runpy import run_path
import sys

from flint import arb, arb_mat


def verify(root, first):
    geometry = run_path(str(first))["geometry"]
    H,Hg,terms,_,edge_groups,groups,fixed,ranges = geometry(root)
    native = run_path(str(root/"ch39-wilson-conductance-law.py"))["data"]
    _,shell,_,edges,base,faces,models = native(root)
    records = loads(Path(__file__).with_suffix(".json").read_text())
    baseline = 2*sum(base.values())
    charges = (2,4,6)
    # 给独立的轻根选取逐支求最小 Q/T²。三条右端由此改写成相对于一阶矩的锥。
    alpha = [min((sum(x*x for x in choice)/sum(choice)**2).lower()
                 for choice in models[q][2]) for q in charges]
    print("三荷目标二阶迹与一阶迹平方之比的共同下端",*[x.str(15) for x in alpha],flush=True)
    pivots = (43,60,61)
    core = arb_mat([[Hg[i,j] for j in pivots] for i in range(3)])
    assert not core.det().contains(0)
    candidates = list(zip(map(str,range(8)),product(*ranges)))+[("selected",fixed)]
    worst=arb(0);floor=arb(1000)
    for key,target in candidates:
        coeff=[arb(n)/10**7 for n in records[key]]
        residual=arb_mat([[v] for v in target])-Hg*arb_mat([[v] for v in coeff])
        for j,correction in zip(pivots,core.solve(residual).entries()):
            coeff[j]+=correction
        face=[arb(0)]*len(terms)
        for orbit,v in zip(groups,coeff):
            for f in orbit: face[f]=v
        assert all(z.contains(0) for z in (H*arb_mat([[v] for v in face])-arb_mat([[v] for v in target])).entries())
        ratios=[[arb(1)]*len(edges) for _ in range(6)]
        for v,(arc,powers) in zip(face,terms):
            for qi,power in enumerate(powers):
                for e in arc: ratios[qi][e]+=v*power
        minimum=min(t.lower() for row in ratios for t in row)
        assert minimum>arb(".20")
        moments=[]
        for k,q in enumerate(charges):
            conduct=[base[e]*v for e,v in zip(edges,ratios[q-1])]
            degrees=[arb(0)]*19
            for (i,j),v in zip(edges,conduct):
                degrees[i]+=v
                degrees[j]+=v
            total=sum(degrees)
            assert (total-baseline-target[k]).contains(0)
            square=sum(v*v for v in degrees)+2*sum(v*v for v in conduct)
            limit=alpha[k]*total**2 if key!="selected" else sum(v*v for v in models[q][3])
            moments.append(square/limit)
        maximum=max(v.upper() for v in moments)
        if key!="selected":
            floor=min(floor,minimum);worst=max(worst,maximum)
        print("角点或已选支、六荷最小边比、三个二阶矩与上限之比",
              key,minimum.str(11),*[v.str(13) for v in moments],flush=True)
        assert maximum < (arb(".982") if key!="selected" else arb(".951"))
    assert floor>arb(".22156") and worst<arb(".981303")
    print("全盒六荷最小边比下界与二阶矩统一比上界",floor.str(16),worst.str(16),flush=True)
    lift_complete(root,records["lift"],shell,edges,base,faces,terms,groups,edge_groups)


def lift_complete(root,record,shell,edges,base,faces,terms,groups,edge_groups):
    towers=run_path(str(root/"ch39-three-even-positive-towers.py"))
    # 原逆谱证书给出三组各有 19 个未知电导的精确解，及其明确盒宽。
    towers["verify"](root)
    labels=run_path(str(root/"ch39-mincut-real-triangle-tomography.py"))["enumerate_faces"]()
    def key(i):
        side=faces[i][0]
        clock,(a,b)=labels[side[0][0],side[0][1],side[1][1]]
        return (clock,a,abs(b) if clock==(0,0) else b)
    ratios=[[{} for _ in edges] for _ in range(6)]
    for j,g in enumerate(groups):
        for f in g:
            es,powers=terms[f]
            for e in es:
                for q,z in enumerate(powers):
                    ratios[q][e][j]=ratios[q][e].get(j,arb(0))+z
    # 对每条实际边比较入射面的整数相位标签多重集，保证代表边方程覆盖整个三边轨道。
    for g in groups:
        for orb in edge_groups:
            tags=[sorted(key(f) for f in g if e in terms[f][0]) for e in orb]
            assert all(tag==tags[0] for tag in tags)
    def sub(charges, indices):
        return arb_mat([[ratios[q-1][orb[0]].get(i,0) for i in indices] for q in charges for orb in edge_groups])
    I=record["pivots"]
    J=sorted(set(range(len(groups)))-set(I))
    assert len(I)==len(set(I))==135 and len(J)==44
    free=dict(record["free"])
    assert set(free)<=set(J)
    b=arb_mat([[arb(free.get(j,0))/10**12] for j in J])
    E=sub((2,4,6),I)
    inverse=E.inv()
    assert all(z.contains(0) for z in (E*inverse-arb_mat([[int(i==j) for j in range(135)] for i in range(135)])).entries())
    T=sub((2,4,6),J)
    O=sub((1,3,5),I)
    P=O*inverse
    original=towers["orbits"](shell["geo"])
    position={e:j for j,g in enumerate(original) for e in g}
    d=[]
    minimum=arb(1)
    for q in (2,4,6):
        weights,active,_=towers["DATA"][q]
        for orb in edge_groups:
            j=position[edges[orb[0]]]
            assert all(position[edges[i]]==j for i in orb)
            pole=shell["soft"]["pole"]; ring=shell["soft"]["star"]["orbit"]
            soft=[pole in edges[i] and any(k in ring for k in edges[i]) for i in orb]
            assert len(set(soft))==1
            v=arb(weights[j])/10**15 + (arb(0,arb("2e-11")) if j in active else 0)
            t=v/base[edges[orb[0]]]
            minimum=min(minimum,t.lower())
            d.append([t-1])
    d=arb_mat(d)
    odd=P*d+(sub((1,3,5),J)-P*T)*b
    lower=min((1+v).lower() for v in odd.entries())
    upper=max((1+v).upper() for v in odd.entries())
    large=max((1+v).lower() for v in odd.entries())
    assert large>arb("1.13e14")
    a=inverse*(d-T*b)
    high=max(abs(v).upper() for v in a.entries()+b.entries())
    print("偶荷原始电导响应的秩、偶荷最小边比、整盒奇荷最小与最大边比、面权最大绝对值上界",
          len(I),minimum.str(12),lower.str(12),upper.str(12),high.str(12),flush=True)
    assert minimum>arb(".125") and lower>arb(".03999") and upper<arb("1.14e14") and high<arb("4.20e14")


if __name__=="__main__":
    root=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(__file__).resolve().parent
    first=Path(sys.argv[2]).resolve() if len(sys.argv)>2 else root/"ch39-signed-wilson-cone.py"
    del sys.argv[1:]
    verify(root,first)
