"""十九点壳上任意轮换不变的非负站点耦合及轻 Higgs 方向的味范数比。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat, ctx


filters=run_path(str(Path(__file__).with_name("ch39-flavor-shell-filters.py")))


def orbit_partition(sites):
    rot=[sites.index((y,z,x)) for x,y,z in sites]
    orbits=[sorted({i,j,rot[j]}) for i,j in enumerate(rot) if i==min(i,j,rot[j])]
    assert len(orbits)==7
    return orbits


def blocks(left,right,H,conjugate):
    orbits=orbit_partition(filters["local"]["local"]["shell"]["geo"]["sites"])
    modes={}
    for _,k,phi in H:
        entries=[]
        for _,a,u in left:
            for _,b,v in right:
                if (a+b+(-k if conjugate else k))%3==0:
                    entries.append([sum(u[i]*v[i]*(phi[i].conjugate() if conjugate else phi[i]) for i in orbit)
                                    for orbit in orbits])
        assert len(entries)==3
        modes[k]=acb_mat(entries)
    return modes


def real_gram(C):
    return arb_mat([[(sum(C[k,i].conjugate()*C[k,j] for k in range(C.nrows()))).real
                     for j in range(C.ncols())] for i in range(C.ncols())])


def stationary_form(M):
    n=M.nrows()
    values=[]
    singular=[]
    for size in range(1,n+1):
        for indices in combinations(range(n),size):
            if size==1:
                values.append((M[indices[0],indices[0]],indices))
                continue
            A=arb_mat([[M[i,j] for j in indices] for i in indices])
            try:
                w=A.solve(arb_mat([[1] for i in indices]))
            except ZeroDivisionError:
                singular.append(indices)
                continue
            norm=sum(w[i,0] for i in range(size))
            assert not norm.contains(0)
            weights=[w[i,0]/norm for i in range(size)]
            if all(s>0 for s in weights):
                values.append((1/norm,indices))
            else:
                assert any(s<0 for s in weights)
    return values,singular


def simple_max(M):
    values, singular=stationary_form(M)
    assert not singular
    return max((v for v,ix in values),key=lambda val:val.upper())


def ratio_bracket(G,K):
    lo,hi=arb(0),arb(16)
    while True:
        value=simple_max(G-hi**2*K)
        if value<0:
            break
        assert value>0
        hi*=2
        assert hi<4096
    for _ in range(23):
        mid=(lo+hi)/2
        value=simple_max(G-mid**2*K)
        if value>0:
            lo=mid
        else:
            assert value<0
            hi=mid
    return lo.union(hi)


def verify():
    reference=run_path(str(Path(__file__).with_name("ch39-shell-flavor-matching-windows.py")))
    ctx.prec=448
    up_ref,down_ref,lepton_ref=[
        sum(y*y for y in reference[name]).sqrt() for name in ("up","down","leptons")
    ]
    scale=arb(".999")/arb("1.001")
    quark_floor=scale*up_ref/down_ref
    lepton_floor=scale*up_ref/lepton_ref
    assert arb("61.36519")<quark_floor<arb("61.36520")
    assert arb("84.98147")<lepton_floor<arb("84.98148")
    print("各汤川中心上下千分之一的上/下与上/带电轻子参考窗",
          quark_floor.str(10), (up_ref/down_ref/scale).str(10),
          lepton_floor.str(10), (up_ref/lepton_ref/scale).str(10))
    for number in (22027,28325):
        h=acb.zeta_zero(number).imag
        shell=filters["local"]["local"]["shell"]
        neutral=shell["laplacian"](0,h)
        center=sum(neutral[i,i].real for i in range(19))/18
        left,up,down,H,lep,charged=[filters["spectra"](q,h,center) for q in (1,-4,2,3,-3,6)]
        U,D=blocks(left,up,H,False),blocks(left,down,H,True)
        E=blocks(lep,charged,H,True)
        sites=shell["geo"]["sites"]
        orbits=orbit_partition(sites)
        singleton=next(i for i,(x,y,z) in enumerate(sites) if x==y==z)
        deficits=[]
        for k in range(3):
            # 非平凡轮换荷的 Higgs 在唯一单点轨道为零，两部门都对此权重不敏感。
            keep=[i for i in range(7) if any(abs(C[k][j,i])>0 for C in (U,D,E) for j in range(3))]
            # 唯一固定点上的非平凡特征向量严格为零，区间值也应含零。
            assert set(range(7))-set(keep)==({orbits.index([singleton])} if k else set())
            assert len(keep)==(7 if k==0 else 6)
            full=[real_gram(C[k]) for C in (U,D,E)]
            G,K,J=[arb_mat([[C[i,j] for j in keep] for i in keep]) for C in full]
            du=ratio_bracket(G,K)
            eu=ratio_bracket(G,J)
            deficits.append(simple_max(G-quark_floor**2*K))
            assert eu<arb("3.34")<lepton_floor
            if number==28325:
                assert du<arb("20.58002")<quark_floor
            print(number,k,"非负轨道耦合的上/下、上/带电轻子范数比上确界",
                  du.str(9),eu.str(9))
        if number==22027:
            assert deficits[0]<0 and deficits[1]>0 and deficits[2]<0
            miss=min(-deficits[0].upper(),-deficits[2].upper())
            leak=deficits[1].upper()/(miss+deficits[1].upper())
            assert leak<arb(".000152")
            print("命中夸克比时，其他两个特征标合计最多份额",leak.str(9))
        else:
            assert all(v<0 for v in deficits)


if __name__=="__main__":
    verify()
