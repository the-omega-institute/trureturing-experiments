"""不标定味系数的局域单 Higgs 折叠，用三站点子式检验三代质量乘积。"""
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


shell = run_path(str(Path(__file__).with_name("ch39-three-ray-chiral-boundary.py")))


def light_basis(q, height, center):
    lap = shell["laplacian"](q,height)
    eig, vec = lap.eig(right=True,algorithm="rump")
    selected = []
    for j,z in enumerate(eig):
        h = abs(z.real*(z.real-center))
        assert h < 9 or h > 9
        if h < 9:
            norm = sum(abs(vec[i,j])**2 for i in range(19)).sqrt()
            selected.append([vec[i,j]/norm for i in range(19)])
    assert len(selected)==3 and all(not z.overlaps(w) for j,z in enumerate(eig) for w in eig[:j])
    return acb_mat(list(zip(*selected)))


def minors(matrix):
    return {s:acb_mat([[matrix[i,j] for j in range(3)] for i in s]).det()
            for s in combinations(range(19),3)}


def verify():
    ctx.prec = 384
    # 与先前 1 TeV 参考相同，六个 Yukawa 各允许千分之一相对改动。
    central = (arb(".00000615")*arb(".00311")*arb(".8616")
               /(arb(".0000135")*arb(".000268")*arb(".01401")))
    low = central*(arb(".999")/arb("1.001"))**3
    high = central*(arb("1.001")/arb(".999"))**3
    print("目标上下型质量乘积比", central.str(9), low.str(9), high.str(9))
    for number in (22027,28325):
        height = acb.zeta_zero(number).imag
        neutral = shell["laplacian"](0,height)
        center = sum(neutral[i,i].real for i in range(19))/18
        up = minors(light_basis(-4,height,center))
        down = minors(light_basis(2,height,center))
        left = minors(light_basis(1,height,center))
        ratios = []
        for s in up:
            assert abs(down[s])>0
            ratios.append(abs(up[s])/abs(down[s]))
        target = low.union(high)
        feasible = [r for r in ratios if r.overlaps(target)]
        assert len(ratios)==969 and not feasible
        cap = max(x.upper() for x in ratios)/low
        print(number, "三站点支持数、最小最大比、能进入目标窗的支持数", len(ratios),
              min(x.lower() for x in ratios).str(9),max(x.upper() for x in ratios).str(9),len(feasible))
        print("任意支持时下型行列式模与其子式项模和之比上限", cap.upper().str(9))
        # 在四个指定站点上，前三个 Higgs 分量取一，末个取 z。
        # 令上型行列式等于中心比乘下型行列式的共轭。
        sites = (0,1,2,3)
        cubic = list(combinations(sites,3))
        a = [left[s]*up[s] for s in cubic]
        b = [(left[s]*down[s]).conjugate() for s in cubic]
        den = sum(a[1:])-central*sum(b[1:])
        assert abs(den)>0
        z = (central*b[0]-a[0])/den
        p, q = a[0]+z*sum(a[1:]), b[0]+z*sum(b[1:])
        assert abs(q)>0 and abs(z)>0 and target.contains_interior(abs(p/q))
        print("四点构造的末站模长、上/下行列式模比", abs(z).str(9),abs(p/q).str(9))


if __name__ == "__main__":
    verify()
