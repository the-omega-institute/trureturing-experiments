"""对指定的连续 ζ 零点窗口逐枚核查 Higgs 真空的顶底层级。"""
import json
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def prepare_certificate(case, rotation, vacuum, positive):
    size = len(rotation)
    assert len(case["first"]) == size
    gram = vacuum["gram"](rotation,case)
    positive(gram)
    v = [arb(z)/10**6 for z in case["first"]]
    local = acb_mat([[gram[i,i]+gram[j,j]-2*gram[i,j]+acb(0,v[i]-v[j])
                      for j in range(size)] for i in range(size)])
    total = arb(0)
    extra = acb_mat(size,size)
    for weight, coordinates in case.get("phase_terms", []):
        assert weight > 0
        coeff = arb(weight)/10**7
        total += coeff
        z = [acb(re,im)/10**5 for re,im in coordinates]
        assert len(z) == size and all(abs(a) <= 1 for a in z)
        for i in range(size):
            for j in range(size):
                if i != j:
                    orig, dest = i,j
                    mean = acb(0)
                    for _ in range(3):
                        mean += z[orig]*z[dest].conjugate()
                        orig,dest = rotation[orig],rotation[dest]
                    extra[i,j] += coeff*mean/3
                else:
                    extra[i,j] += coeff
    return local,total,extra


def lowest_projector(matrix,bases):
    proj = acb_mat(19,19)
    for basis in bases:
        block = basis.conjugate().transpose()*matrix*basis
        values,vectors = block.eig(right=True)
        below = [j for j, value in enumerate(values) if value.real < 9]
        assert len(below) == 1 and all(value.real > 9 or value.real < 9 for value in values)
        idx, = below
        v = basis*acb_mat([[vectors[i,idx]] for i in range(len(values))])
        vv = v*v.conjugate().transpose()
        proj += vv/sum(abs(row[0])**2 for row in v.tolist())
    return proj


def verify():
    root = Path(__file__).parent
    ctx.prec = 192
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    vacuum = run_path(str(root/"ch39-higgs-phase-vacuum.py"))
    positive = run_path(str(root/"ch39-three-mode-higgs-hierarchy.py"))["positive"]
    common = run_path(str(root/"ch39-all-common-shell-ranks.py"))
    ctx.prec = 192
    sites = shell["geo"]["sites"]
    rotation = [sites.index((y,z,x)) for x,y,z in sites]
    data = json.loads((root/"ch39-zero-window-vacuum.json").read_text())
    # 导出器把有理 Gram 与一阶系数统一为分母 10^7。
    certs = []
    for entry in data["certificates"]:
        case = dict(entry,gram=[arb(n)/10 for n in entry["gram"]],
                    first=[arb(n)/10 for n in entry["first"]])
        certs.append(prepare_certificate(case,rotation,vacuum,positive))
    bases = list(common["character_bases"](shell))
    assert len(data["choices"]) == len(data["multipliers"]) == 1238
    heights = acb.zeta_zeros(10855,1238)
    assert len(heights) == 1238 and all(a.imag < b.imag for a,b in zip(heights,heights[1:]))
    local_count = 0
    for zero,choice,coef in zip(heights,data["choices"],data["multipliers"]):
        assert 0 <= choice < len(certs) and coef > 0
        local,total,extra = certs[choice]
        height = zero.imag
        lap = shell["laplacian"](3,height)
        p1,p4,p2 = [lowest_projector(shell["laplacian"](q,height),bases) for q in (1,4,-2)]
        matrix = acb_mat(19,19)
        for i in range(19):
            for j in range(19):
                down = p1[i,j]*p2[i,j].conjugate()
                up = p1[i,j].conjugate()*p4[i,j]
                matrix[i,j] = (local[i,j]+total)*lap[i,j]-extra[i,j]*lap[i,j].conjugate()
                matrix[i,j] += coef*(down-up/300)
        positive(matrix)
        local_count += int(total.is_zero())
    print("窗口、局部极小证书数、使用全局能量比较的编号数、统一顶底上界",
          (10855,12092),local_count,1238-local_count,30)


if __name__ == "__main__":
    verify()
