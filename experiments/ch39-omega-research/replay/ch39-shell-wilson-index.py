"""双钟球壳的三边界、Wilson 质量窗口与小方格规范曲率界。"""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
from runpy import run_path


def verify_clock_boundary(shell):
    from collections import Counter
    edges = shell["geo"]["edges"]
    clock_edges = [(i,j) for i,j,d,_ in edges if d in (724,804)]
    incidence = Counter(v for edge in clock_edges for v in edge)
    inside = set(shell["boundary_sites"]())
    assert len(clock_edges) == 12 and len(inside) == 3
    assert {v for v in incidence if incidence[v] == 4} == inside
    assert max(incidence.values()) == 4 and all(len(inside.intersection(e)) == 1 for e in clock_edges)
    inside_sites = {tuple(shell["geo"]["sites"][v]) for v in inside}
    assert inside_sites == {(23,13,13), (13,23,13), (13,13,23)}
    for v in inside:
        neighbors = {j if i == v else i for i,j,_,_ in edges if v in (i,j)}
        assert len(neighbors) == 18 and len(neighbors.intersection(inside)) == 2
    print("时钟边、唯一三点覆盖及内块 Gershgorin 下界", len(clock_edges), sorted(inside), 18-2)


def verify_mass_window():
    low, high = F(8685,1000), F(96,10)
    center, scale = (low+high)/2, high/2
    gap = (high-low)/high
    mass = lambda level: (level-center)/scale
    assert gap == F(61,640)
    assert mass(0) == -2+gap and mass(low) == -gap and mass(high) == gap
    # 四维自由符号中相异方向的二次交叉项保留为多仿射式。
    sample = [F(1,5), F(1,2), F(2,3), F(9,5)]
    for m in (mass(0), mass(low), F(-1), mass(high), mass(36)):
        direct = sum(u*(2-u) for u in sample) + (m+sum(sample))**2
        expanded = m*m+2*(1+m)*sum(sample)+2*sum(a*b for a,b in combinations(sample,2))
        assert direct == expanded
        assert min(abs(m+2*k) for k in range(5)) >= gap
    # 每对方向的三个误差界分别为 2、1、2 倍的面范数。
    loss = 6*(2+1+2)
    limit = (1-(1-gap)**2)/loss
    assert limit == F(74359,12288000)
    assert (1-gap+F(17,1000))**2 < 1-loss*F(1,200)
    print("自由隙、面曲率上限、面误差不超过 1/200 时的隙下界", gap, limit, ">0.017")


def verify_flux_example():
    size = 96
    modulus = size**2

    def links(x,y):
        return (-size*y if x == size-1 else 0), x

    for x in range(size):
        for y in range(size):
            h,v = links(x,y)
            next_v = links((x+1)%size,y)[1]
            next_h = links(x,(y+1)%size)[0]
            assert (h+next_v-next_h-v)%modulus == 1
    # 两个不相交坐标平面拷贝此场，混合面的链相互抵消。
    assert 2*F(22,7)*6/modulus < F(1,200)
    print("周期边长、每个活跃平面的单位通量", size, 1)


if __name__ == "__main__":
    root = Path(__file__).parent
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    verify_clock_boundary(shell)
    run_path(str(root/"ch39-all-zero-three-mode-gap.py"))["verify"]()
    verify_mass_window()
    verify_flux_example()
