"""验证十九点壳的手征表示、两枚真实零点的换点余量及谱线计数。"""
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


def check_representation():
    left_doublets = [1, 1, 1, -3]
    right_singlets = [4, 4, 4, -2, -2, -2, -6]
    left_hypercharges = [y for y in left_doublets for _ in range(2)]
    conjugated = left_hypercharges + [-y for y in right_singlets]
    assert sum(left_doublets) == 0 and len(left_doublets) % 2 == 0
    assert len(conjugated) == 15 and sum(conjugated) == sum(q**3 for q in conjugated) == 0
    assert Counter(abs(q) for q in conjugated if q%2) == {1:6, 3:2}
    # 外部独立规范场的一个充分界。对双重态用乘积方格的三角不等式。
    assert max(F(1,400) + abs(q)*F(1,2400) for q in left_doublets) < F(1,200)
    assert max(abs(q)*F(1,2400) for q in right_singlets) < F(1,200)


def check_switching():
    root = Path(__file__).parent
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    geo = shell["geo"]
    for d in (724, 804):
        pairs = [(i,j) for i,j,edge_d,_ in geo["edges"] if edge_d == d]
        assert Counter(i for i,_ in pairs) == {i:1 for i in range(6)}
        assert Counter(j for _,j in pairs) == {6:2, 7:2, 8:2}

    ctx.prec = 192
    center, scale = arb("9.1425"), arb("4.8")
    heights = [acb.zeta_zero(n).imag for n in (11431,133931)]
    clocks = [geo["phases"](h) for h in heights]
    first, second = clocks
    chord = sum(abs(first[d]-second[d]) for d in (724,804))
    assert arb(".00976084") < chord < arb(".00976085")
    curvature_loss = 1-(1-arb(30)/200).sqrt()
    assert arb(".1049") - curvature_loss - (6*arb(2).sqrt()/100 + 6*arb("1e-8"))/scale > arb(".009")
    for q in (1,2,3,4,6):
        lap = [shell["laplacian"](q,h) for h in heights]
        masses = sorted(((z.real-center)/scale for z in lap[0].eig()),key=float)
        free_gap = min(min(abs(m+2*k) for k in range(5)) for m in masses)
        assert free_gap > arb(".1049")
        assert sum(-2 < m < 0 for m in masses) == 3 and sum(m > 0 for m in masses) == 16
        bound = (arb(2).sqrt()*sum(abs(first[d]**q-second[d]**q) for d in (724,804))
                 + 6*arb("1e-8"))/scale
        assert bound < arb(".017254") and free_gap-curvature_loss-bound > arb(".0096")
        # 周期格为 2×1×1×1 时，平凡外场的平移自伴，所有反对称跳项为零。
        # 两格交替放置两种壳质量，Wilson 实部给站点对角一、跨点负一。
        alternating = acb_mat(38,38)
        for site, block in enumerate(lap):
            for i in range(19):
                n = 19*site+i
                for j in range(19):
                    alternating[n,19*site+j] = block[i,j]/scale
                alternating[n,n] += 1-center/scale
                alternating[n,19*(1-site)+i] = -1
        values = [z.real for z in alternating.eig()]
        assert sum(v < 0 for v in values) == 3 and sum(v > 0 for v in values) == 35
        assert min(abs(v) for v in values) > arb(".1")
        print("荷、参考自由隙、换点质量算子差上界",q,free_gap.str(12),bound.str(12))
    print("共同两钟弦距之和",chord.str(15))


if __name__ == "__main__":
    check_representation()
    check_switching()
