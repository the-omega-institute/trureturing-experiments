"""三套实际重层质量独立受扰时，共同非负累计引力轮廓的全区间分离。"""
from pathlib import Path
from runpy import run_path

from flint import arb


def verify(root):
    _, fields=run_path(str(root/"ch39-three-index-free-atoms.py"))["spectra"](root)
    dual=run_path(str(root/"ch39-three-index-positive-profile-obstruction.py"))["coefficients"]
    api=run_path(str(root/"ch39-common-cutoff-negative-mass.py"))
    rho=(arb("246.22")/arb("1.22089e19"))**2/2
    limit=arb("4e18")
    displacement=rho*limit
    edge,cut=arb("1e-14"),arb("1e-5")
    # 每组最后 32 行是 16 对同质量起点的 Dirac 伴侣与标量。
    per_field=[]
    for field in fields:
        heavy=list(field[4:])
        assert len(heavy)==114
        for i in range(82,114,2):
            g,a,z=heavy[i]
            g1,a1,m=heavy[i+1]
            assert g==-8 and g1==4 and (a-a1).contains(0) and abs(m-z)<rho
            heavy[i+1]=(g1,a1,z)
        per_field.append(field[:4]+heavy)
    def moving(i):
        return i>=83 and i%2==1
    reference=per_field[2][4+82][2]
    def interval_mass(i,z):
        return z+arb(0,displacement.upper()) if moving(i) else z
    pairs=[]
    for field in per_field:
        bare=field[4+82][2]
        low=bare.union(bare-displacement)
        assert low>0
        assert all(interval_mass(i,z)>bare for i,(_,_,z) in enumerate(field[4:]) if i not in (82,83))
        pairs.append(low)
    assert pairs[0]-pairs[2]>arb("1.747e-13") and pairs[1]-pairs[0]>arb("7.7e-14")
    print("N-1 减中心、中心减 N+1 的最轻重行质量平方比间隔下界",
          (pairs[1]-pairs[0]).lower().str(14),(pairs[0]-pairs[2]).lower().str(14),flush=True)
    normal=list(map(arb,("1","1",".47",".495")))
    radius=arb(".06")
    directions=[[v+(sign*radius if i==j else 0) for i,v in enumerate(normal)]
                for j in range(4) for sign in (-1,1)]
    tails=[]
    counts=[]
    for direction in directions:
        coeff=dual(per_field,direction)
        positive=[]
        tail=arb(0)
        for n,i,z,k,w in coeff:
            if n==2 and i in (82,83):
                positive.append((k+arb(0).union(edge)*w).lower())
            else:
                # 实际最轻行至多为它的固定伴侣质量，因此可用 reference 算保守间隙。
                gap=interval_mass(i,z)-reference
                assert gap>arb("1.747e-13")
                tail+=(abs(k)+edge*abs(w))*(-gap/edge).exp()
        tails.append(min(positive)-tail)
        assert tails[-1]>arb(".27")

        # 裸谱中同一质量的两行与抽出的参考指数完全相消。
        rows=[(arb(0) if n==2 and i in (82,83) else z-reference,k,w) for n,i,z,k,w in coeff]
        scalars=[(z-reference,k,w) for _,i,z,k,w in coeff if moving(i)]
        assert len(scalars)==48
        large=api["long_series"](rows,27)
        def error(x):
            # 固定抽出 exp(-reference/x)。均值定理作用在完整质量扰动路径。
            return displacement*sum(abs(k/x+w)*((displacement-d)/x).exp() for d,k,w in scalars)
        count=0
        for l,r,fn in ((edge,cut,lambda x:api["small_derivative"](x,rows)),(cut,arb(1),large)):
            todo=[(l,r)]
            while todo:
                left,right=todo.pop()
                box=left.union(right)
                if fn(box)>error(box):
                    count+=1
                    continue
                assert right-left>arb("1e-24"),(direction,left,right,fn(box),error(box))
                mid=((left+right)/2).mid()
                todo.extend(((left,mid),(mid,right)))
        counts.append(count)
    print("门户绝对幅度、质量平方比的最大绝对改变量、零端最小分离下界",
          limit.str(10),displacement.upper().str(14),min(tails).lower().str(14),flush=True)
    print("八个分离方向的连续覆盖闭段数",counts,flush=True)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
