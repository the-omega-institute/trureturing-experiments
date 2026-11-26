"""同一实际壳谱与一跳支撑下，荷一到荷二可以实现任意三阶酉味映射。"""
from itertools import permutations
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def verify():
    root = Path(__file__).resolve().parent
    light = run_path(str(root / "ch39-local-yukawa-plucker.py"))
    prev = run_path(str(root / "ch39-joint-zero-diagonal-light-filters.py"))
    adjoint = run_path(str(root / "ch39-joint-zero-shell-local-transfer.py"))["adjoint"]
    ctx.prec = 320
    height = acb.zeta_zero(2_306_888_100_085).imag
    shell = light["shell"]
    delta = arb(23)/25*arb(2)**-20/(height**2+shell["soft"]["point"]["study"]["offset"])
    center = 15-delta/3
    u = light["light_basis"](1,height,center)
    v = light["light_basis"](2,height,center)
    sites = range(19)
    identity = acb_mat([[int(i==j) for j in sites] for i in sites])
    p = identity-u*adjoint(u)
    q = identity-v*adjoint(v)
    edges = {(i,i) for i in sites}
    for i,j,_,_ in shell["geo"]["edges"]:
        edges.update(((i,j),(j,i)))
    blocked = [(i,j) for i in sites for j in sites if (i,j) not in edges]
    assert len(blocked)==72
    h = acb_mat([[q[i,k]*p[l,j] for k,l in blocked] for i,j in blocked])
    prev["certify_positive"](h)
    pairs = [(a,b) for a in range(3) for b in range(3)]
    design = acb_mat([[v[i,a]*u[j,b].conjugate() for a,b in pairs] for i,j in blocked])
    gram = adjoint(design)*h.solve(design)
    # 对任意复 M，修正块的 Frobenius 范数平方等于 vec(M)*gram vec(M)。
    # 逐对应用 2|z_i z_j| <= |z_i|²+|z_j|²，得到对角二次型上界。
    row_bound = {}
    for i, ab in enumerate(pairs):
        row_bound[ab] = gram[i,i].real.upper()+sum(abs(gram[i,j]).upper() for j in range(9) if i!=j)
    # 酉矩阵条目模方组成双随机矩阵。线性泛函最大值在六个置换之一取得。
    bounds = [sum(row_bound[a,b] for a,b in enumerate(perm)) for perm in permutations(range(3))]
    worst = max(v.upper() for v in bounds)
    print("六种置换上界",*[z.str(12) for z in bounds])
    print("对任意酉味矩阵的补块平方范数上界",worst.str(16))
    assert worst < arb(".986")


if __name__ == "__main__":
    verify()
