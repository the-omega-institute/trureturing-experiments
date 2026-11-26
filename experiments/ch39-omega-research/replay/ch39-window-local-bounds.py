"""逐枚验证真实零点窗口内局部 Higgs 真空的顶底比和相位变分界。"""
import json
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def column(entries):
    assert len(entries) == 19
    return acb_mat([[acb(real,imag)/10**7] for real,imag in entries])


def square(vector):
    return vector * vector.conjugate().transpose()


def phase_test(lap, flavor, gram, first, cut, positive):
    # 单位方向的梯度及 Hessian 以 19 个实站点相角的欧氏范数计。
    assert len(first) == 19
    first = [arb(a)/10**7 for a in first]
    positive(gram)
    allowance = (sum(gram[i,i].real for i in range(19)) + sum(abs(a) for a in first))/10**9
    matrix = acb_mat(flavor)
    for i in range(19):
        for j in range(19):
            factor = gram[i,i] + gram[j,j] - 2*gram[i,j] + acb(0,first[i]-first[j])
            matrix[i,j] += factor*lap[i,j]
        matrix[i,i] -= allowance
    positive(matrix + cut)


def verify_cover(data, lap, flavor, positive):
    a, b = map(column,data["axes"])
    total = square(a) + square(b)
    def visit(node, i, m, j, s):
        if isinstance(node,list):
            assert len(node) == 4
            return sum(visit(node[2*u+v],2*i+u,2*m,2*j+v,2*s)
                       for u in range(2) for v in range(2))
        theta = (2*i+1)*arb.pi()/(4*m)
        phi = (2*j+1)*arb.pi()/s
        epsilon = arb(node["epsilon"])/10**7
        assert (arb.pi()/(4*m)).sin()**2 + (arb.pi()/(2*s)).sin()**2 <= epsilon
        coefficient = arb(node["tau"])/10**7
        assert coefficient >= 0
        normal = a*(-theta.sin()) + b*(acb(0,phi).exp()*theta.cos())
        cut = (square(normal)-total*epsilon)*coefficient
        gram = acb_mat(19,19)
        entries = iter(node["gram"])
        for u in range(19):
            for v in range(u,19):
                gram[u,v] = gram[v,u] = arb(next(entries))/10**7
        assert next(entries,None) is None
        phase_test(lap,flavor,gram,node["first"],cut,positive)
        return 1
    assert len(data["roots"]) == 18
    return sum(visit(data["roots"][6*i+j],i,3,j,6) for i in range(3) for j in range(6))


def verify():
    root = Path(__file__).parent
    ctx.prec = 192
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    window = run_path(str(root/"ch39-window-phase-vacua.py"))
    vacuum = run_path(str(root/"ch39-higgs-phase-vacuum.py"))
    positive = run_path(str(root/"ch39-three-mode-higgs-hierarchy.py"))["positive"]
    common = run_path(str(root/"ch39-all-common-shell-ranks.py"))
    ctx.prec = 192
    bases = list(common["character_bases"](shell))
    sites = shell["geo"]["sites"]
    rotation = [sites.index((y,z,x)) for x,y,z in sites]
    old = json.loads((root/"ch39-zero-window-vacuum.json").read_text())
    cases = json.loads((root/"ch39-window-local-bounds.json").read_text())
    heights = acb.zeta_zeros(10855,1238)
    assert len(heights) == len(old["choices"]) == len(old["multipliers"]) == 1238
    assert all(a.imag < b.imag for a,b in zip(heights,heights[1:]))
    used, covers, leaves = set(), [], 0
    for n,zero in enumerate(heights,10855):
        i = n-10855
        previous = old["certificates"][old["choices"][i]]
        if previous.get("phase_terms"):
            case, weight = cases[str(n)], 300
            used.add(str(n))
        else:
            case, weight = previous, old["multipliers"][i]
        assert weight > 0
        lap = shell["laplacian"](3,zero.imag)
        p1, p4, pm2 = [window["lowest_projector"](shell["laplacian"](q,zero.imag),bases)
                       for q in (1,4,-2)]
        flavor = acb_mat([[weight*(p1[i,j]*pm2[i,j].conjugate()-p1[i,j].conjugate()*p4[i,j]/300)
                           for j in range(19)] for i in range(19)])
        if "roots" in case:
            covers.append(n)
            leaves += verify_cover(case,lap,flavor,positive)
        else:
            gram = vacuum["gram"](rotation,{"gram":[arb(v)/10 for v in case["gram"]]})
            phase_test(lap,flavor,gram,case["first"],acb_mat(19,19),positive)
    assert used == set(cases)
    print("零点总数、单份局部证书数、分区证明的编号及总叶数、顶底比上界、变分尺度",
          len(heights),len(heights)-len(covers),covers,leaves,30,"1e-9")


if __name__ == "__main__":
    verify()
