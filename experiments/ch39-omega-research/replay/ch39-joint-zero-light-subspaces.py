"""联合近真空零点的五荷轻空间、逐站相位可实现的重叠及其规范无关下界。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


def monomial_bound(p, r):
    # 同一置换作用在行与列。先放开每行的列置换，再对整行的匹配作认证上界。
    rows = []
    for mat in (p, r):
        rows.append([sorted((abs(mat[i,j]).abs_upper() for j in range(19)), reverse=True)
                     for i in range(19)])
    # 任意置换的各行内积至多为排序后的内积。
    sizes = [[sum(a*b for a,b in zip(u,v,strict=True)) for v in rows[1]] for u in rows[0]]
    # 用整数向上取整，求指派问题的对偶可行势。每个比较都可在整数上检查。
    scale = 10**12
    cost = [[-int((v*scale).upper().ceil().unique_fmpz()) for v in row] for row in sizes]
    n = len(cost)
    row_potential, col_potential = [0]*(n+1), [0]*(n+1)
    match, previous = [0]*(n+1), [0]*(n+1)
    for row in range(1,n+1):
        match[0] = row
        column = 0
        best = [None]*(n+1)
        used = {0}
        while match[column]:
            i = match[column]
            for j in range(1,n+1):
                if j in used:
                    continue
                reduced = cost[i-1][j-1]-row_potential[i]-col_potential[j]
                if best[j] is None or reduced<best[j]:
                    best[j],previous[j] = reduced,column
            next_col = min((j for j in range(1,n+1) if j not in used),key=lambda j:best[j])
            delta = best[next_col]
            for j in range(n+1):
                if j in used:
                    row_potential[match[j]] += delta
                    col_potential[j] -= delta
                else:
                    best[j] -= delta
            column = next_col
            used.add(column)
        while column:
            j = previous[column]
            match[column] = match[j]
            column = j
    assert all(row_potential[i+1]+col_potential[j+1]<=v for i,row in enumerate(cost) for j,v in enumerate(row))
    upper = -sum(row_potential[1:])-sum(col_potential[1:])
    return 1-arb(upper)/(3*scale)


def verify():
    root = Path(__file__).resolve().parent
    light = run_path(str(root / "ch39-local-yukawa-plucker.py"))
    shell = light["shell"]
    ctx.prec = 320
    h = acb.zeta_zero(2_306_888_100_085).imag
    delta = arb(23)/25*arb(2)**-20/(h*h+shell["soft"]["point"]["study"]["offset"])
    center = 15-delta/3
    charges = (1, 2, -4, -3, 6)
    bases = {q:light["light_basis"](q, h, center) for q in charges}
    density = {q:[sum(abs(b[i,j])**2 for j in range(3))/3 for i in range(19)]
               for q,b in bases.items()}
    projectors = {q:frame*frame.conjugate().transpose() for q,frame in bases.items()}
    for q in charges:
        lap = shell["laplacian"](q, h)
        cutoff = 9+(1-delta)**2
        masses = [abs(z.real*(z.real-center)) for z in lap.eig(multiple=True)]
        assert sum(v<9 for v in masses) == 3 and all(v<9 or v>cutoff for v in masses)
    for q, r in combinations(charges, 2):
        overlap = bases[q].conjugate().transpose()*bases[r]
        squared = overlap*overlap.conjugate().transpose()
        eigen = sorted((v.real for v in squared.eig(algorithm="rump")), key=float)
        leakage = 1-sum(eigen)/3
        variation = sum(abs(a-b) for a,b in zip(density[q],density[r],strict=True))/2
        absolute_sum = sum(abs(projectors[q][i,j]*projectors[r][j,i]) for i in range(19) for j in range(19))
        phase_bound = 1-absolute_sum/3
        print("两荷、固定站点基平均泄漏、站点概率和逐元迹给出的下界",
              q, r, leakage.str(11), (variation**2).str(11), phase_bound.str(11))
        assert min(eigen)>0 and max(eigen)<1 and variation>0
        assert phase_bound > arb(".00053")
        monomial = monomial_bound(projectors[q], projectors[r])
        print("允许站点置换与任意相位的下界",monomial.str(11))
        assert monomial>arb(".00042")
        if (q,r)==(1,2):
            assert monomial>arb(".31060011")
        if (q,r)==(1,-4):
            assert monomial>arb(".31097727")
        if (q,r)==(1,2):
            assert leakage>arb(".51702520") and phase_bound>arb(".38321152")
        if (q,r)==(1,-4):
            assert leakage>arb(".52481329") and phase_bound>arb(".37290383")
            print("上型共轭与左手夸克的主角余弦平方",*[v.str(11) for v in eigen])


if __name__ == "__main__":
    verify()
