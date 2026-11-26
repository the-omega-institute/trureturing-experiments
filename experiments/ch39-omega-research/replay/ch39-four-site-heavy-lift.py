"""用保留三个轻波函数的四站点余子式支撑临界重质量提升。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def light_minors(light, triangles):
    return {t:acb_mat([[light[i,j] for j in range(3)] for i in t]).det() for t in triangles}


def cofactor_row(cell, minors):
    return [(-1)**k*minors[cell[:k]+cell[k+1:]] for k in range(4)]


def check_exterior_changes(lap, charge, cells, triangles, original, matching):
    phase = acb(9999,200)/10001
    separation = {cell:arb(0) for cell in cells}
    for i,j in matching:
        perturbed = acb_mat(lap)
        perturbed[i,j] = lap[i,j]*phase**charge
        perturbed[j,i] = perturbed[i,j].conjugate()
        values,vectors = perturbed.eig(right=True)
        order = sorted(range(19),key=lambda k:float(values[k].real))
        assert values[order[2]].real < 9 < values[order[3]].real
        light = acb_mat([[vectors[k,c] for c in order[:3]] for k in range(19)])
        minors = light_minors(light,triangles)
        for cell in cells:
            if i in cell or j in cell:
                continue
            a,b = cofactor_row(cell,original),cofactor_row(cell,minors)
            wedge = [a[u]*b[v]-a[v]*b[u] for u,v in combinations(range(4),2)]
            squared = lambda x: sum((z*z.conjugate()).real for z in x)
            value = squared(wedge)/(squared(a)*squared(b))
            separation[cell] = separation[cell].max(value.lower())
    angle = min(separation.values()).sqrt()
    assert angle > 0
    return angle


def verify():
    directory = Path(__file__).parent
    shell = run_path(str(directory/"ch39-three-ray-chiral-boundary.py"))
    algebra = run_path(str(directory/"ch39-three-mode-higgs-hierarchy.py"))
    adjoint, positive, norm = (algebra[k] for k in ("adjoint","positive","norm"))
    ctx.prec = 320
    height = acb.zeta_zero(11431).imag
    edges = {(i,j) for i,j,_,_ in shell["geo"]["edges"]}
    triangles = list(combinations(range(19),3))
    cells = [s for s in combinations(range(19),4) if all(e in edges for e in combinations(s,2))]
    assert len(triangles) == 969 and len(cells) == 1327
    matching = [(0,2),(1,3),(4,5),(6,7),(8,9)]
    assert all(pair in edges for pair in matching)
    assert len({vertex for pair in matching for vertex in pair}) == 10
    zero_levels, lifted_levels = [], []
    for q in (1,2,3,4,6):
        lap = shell["laplacian"](q,height)
        levels, vectors = lap.eig(right=True)
        order = sorted(range(19),key=lambda j:float(levels[j].real))
        e = [levels[j].real for j in order]
        assert all(a < b for a,b in zip(e,e[1:]))
        zero_levels.extend(e if q in (1,3) else e[:3])
        if q in (1,3):
            continue
        lifted_levels.extend(e[3:])
        basis = []
        for j in order:
            col = acb_mat([[vectors[i,j]] for i in range(19)])
            basis.append(col/norm(col))
        light = acb_mat([[v[i,0] for v in basis[:3]] for i in range(19)])
        heavy = acb_mat([[v[i,0] for v in basis[3:]] for i in range(19)])
        minors = light_minors(light,triangles)
        assert min(abs(z).lower() for z in minors.values()) > arb(".000028")
        lifted = acb_mat(19,19)
        for cell in cells:
            row = cofactor_row(cell,minors)
            for i,a in zip(cell,row):
                for j,b in zip(cell,row):
                    lifted[i,j] += a.conjugate()*b
        assert all(z.contains(0) for row in (lifted*light).tolist() for z in row)
        block = adjoint(heavy)*lifted*heavy
        positive(block)
        lo,hi = {2:(46000,48000),4:(19500,20100),6:(163000,165000)}[q]
        def shifted(strength,ratio):
            result = strength*block
            for j in range(16):
                result[j,j] += e[j+3]-ratio**2
            return result
        positive(shifted(hi,arb("9.661")))
        assert min(v.real for v in shifted(lo,arb("9.659")).eig()) < 0
        print("荷、共同四点权的临界强度包络",q,lo,hi)
        sine = check_exterior_changes(lap,q,cells,triangles,minors,matching)
        assert sine > arb({2:".0000639",4:".0001426",6:".0005905"}[q])
        print("每个四点组均可由一条不接触它的边改变其余子式方向；正弦统一下界",
              sine.str(12))
    zero_levels.sort(key=float)
    lifted_levels.sort(key=float)
    merged = sorted(zero_levels+lifted_levels,key=float)
    assert all(a < b for a,b in zip(merged,merged[1:]))
    gaps = [sum(a < z < b for z in zero_levels) for a,b in zip(lifted_levels,lifted_levels[1:])]
    assert len(zero_levels) == 47 and len(lifted_levels) == 48
    assert sum(n%2 for n in gaps) == 10
    closest = min(abs(z-v).lower() for z in zero_levels for v in lifted_levels)
    assert closest > arb(".00221")
    print("共同谱多项式在实际谱上的最低次数、全正轴非负的最低次数",47+10,2*47)


if __name__ == "__main__":
    verify()
