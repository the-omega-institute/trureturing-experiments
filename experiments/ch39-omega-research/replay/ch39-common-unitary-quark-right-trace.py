"""同一十九维酉输运在右外壳给出的两荷迹下界及固定物理空间的共同弱基覆盖。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def weighted(basis, weights):
    return acb_mat([[basis[i, j]*weights[j] for j in range(basis.ncols())]
                    for i in range(basis.nrows())])*basis.conjugate().transpose()


def eigs(a):
    return sorted((z.real for z in a.eig(algorithm="rump")), key=lambda z: z.mid())


def positive(a):
    for n in range(1, a.nrows()+1):
        principal = acb_mat([[a[i,j] for j in range(n)] for i in range(n)])
        assert principal.det().real > 0


def lower_trace(root):
    common = run_path(str(root/"ch39-common-embedding-central-right-ckm.py"))
    source = 2_306_888_100_086
    _, values, _ = common["verify_nodes"](root, ((source, ".0294"),), ".1", (),
                                          ((6, 7, 10, 13), (6, 7, 10, 12)))
    adj = common["adj"]
    spectrum = run_path(str(root/"ch39-adjacent-pvm-transport.py"))["spectrum"]
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    ctx.prec = 384
    groups = ((0, 3, 4), (1, 2, 5), (6, 7, 8), (9, 15, 16), (10, 13, 17),
              (11, 12, 18), (14,))
    T = acb_mat([[1/arb(len(g)).sqrt() if i in g else 0 for g in groups] for i in range(19)])
    source_ops = []
    measured = []
    char = []
    for charge in (-4, 2):
        _, previous = spectrum(shell, charge, source)
        _, current = spectrum(shell, charge, source+1)
        w = values[source, charge]
        source_ops.append(weighted(previous, w))
        char.append((min(z.lower() for z in w), max(z.upper() for z in w)))
        tr = [abs(current[14, j])**2+abs(sum(current[i, j] for i in (8, 7, 6)))**2/3 for j in range(19)]
        measured.append(adj(T)*weighted(current, tr)*T)
    # 三循环将 H 分成 7+6+6 阶。后两块用各自最小 Gershgorin 左端点。
    h = source_ops[0] + 15*source_ops[1]
    lower = [x.lower() for x in eigs(adj(T)*h*T)]
    sites = shell["geo"]["sites"]
    three = [group for group in groups if len(group) == 3]
    for orient in (-1, 1):
        omega = acb(-arb(1)/2, orient*arb(3).sqrt()/2)
        char_basis = acb_mat(19, 6)
        for j, group in enumerate(three):
            site = group[0]
            for power in range(3):
                char_basis[site,j] = omega**power/arb(3).sqrt()
                x,y,z = sites[site]
                site = sites.index((y,z,x))
            assert site == group[0]
        block = adj(char_basis)*h*char_basis
        edge = min((block[i,i].real - sum(abs(block[i,j]) for j in range(6) if j != i)).lower()
                   for i in range(6))
        lower += [edge]*6
    full = sorted(lower)
    for matrix in measured:
        positive(matrix)
    b = (measured[0].inv() + arb(15)/25*measured[1].inv()).inv()
    vals = eigs(b)
    assert len(vals) == 7 and vals[0]>0 and full[0]>0
    bound = sum(v*w.lower() for v, w in zip(full[:7], reversed(vals), strict=True))
    assert bound>arb(".0125050")
    print("共同酉下的平凡块两迹加权下界", bound.str(15), flush=True)
    return char


def weak_cover(root, char):
    witness = run_path(str(root/"ch39-independent-unitary-quark-right-holdout.py"))
    old = run_path(str(root/"ch39-neighbour-all-frame-support.py"))
    pencil = run_path(str(root/"ch39-triple-geometry-pencil.py"))
    u, v = map(arb, witness["GEOMETRY"][:2])
    _, geom, V, mass, A, B, _ = pencil["setup"](root)
    def test(y, chart, t):
        data = pencil["physical_blocks"](y,chart,t,V,mass,geom,old,())
        if data is None:
            return None
        margins = []
        for sign in (-1,1):
            traces = []
            violations = []
            for (p,s,z,k,w,det),(lo,hi) in zip(data[0],char,strict=True):
                for a in (s-sign*w/u,s+sign*w/v):
                    violations.extend((lo-a, a-hi))
                traces.append(p+(z-u*v*s-(u-v)*sign*w)/(1-u*v))
            sep = max(x.lower() for x in [*violations, arb(".00995")-traces[0]-25*traces[1]])
            if not sep > 0:
                return None
            margins.append(sep)
        return min(margins)
    for chart,(base,other) in enumerate(((A,B),(B,A))):
        stack=[(arb(-1),arb(1),[0,1,2],0)]
        boxes=deep=0
        least=arb(1)
        while stack:
            lo,hi,indices,level=stack.pop()
            boxes+=1
            deep=max(deep,level)
            lines=old["lines"](base,other,lo,hi,indices,geom)
            keep=[]
            for j in indices:
                ok=test(lines[j],chart,lo.union(hi)) if j in lines else None
                if ok is None:keep.append(j)
                else:least=min(least,ok)
            if keep:
                assert level<26,(chart,keep,lo,hi)
                mid=(lo+hi)/2
                stack.extend(((lo,mid,keep,level+1),(mid,hi,keep,level+1)))
        print("铅笔图、检查区间、深度、最低严格分离",chart,boxes,deep,least.str(12),flush=True)


if __name__=="__main__":
    root=Path(__file__).resolve().parent
    ranges = lower_trace(root)
    weak_cover(root, ranges)
    # tr T_u=2；对 H_u 减去谱区间中点后，用算符范数控制两次酉共轭的差。
    gap = (arb(".0125050")-arb(".00995"))/(2*(ranges[0][1]-ranges[0][0]))
    assert gap > arb(".001348")
    print("实现目标时两荷输运的最小射影算符范数距离大于", gap.lower().str(13))
