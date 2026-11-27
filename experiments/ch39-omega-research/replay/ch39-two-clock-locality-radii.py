"""由十九点壳的闭路与顶点重定相验证单位荷两钟联络的两个局域阈值。"""
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb, arb


def check_graph(shell):
    sites = shell["sites"]
    edges = shell["edges"]
    assert len(sites) == 19 and len(edges) == 135
    shortest = min(d for _,_,d,_ in edges)
    assert shortest == 12 and sum(d == shortest for _,_,d,_ in edges) == 3
    assert 1-Fraction(shortest,4*shell["R"]) == Fraction(288,289)
    directed = {}
    inner = set(shell["inner"])
    for i,j,d,_ in edges:
        a,b = int(d == 724),int(d == 804)
        diff = int(sites[j] in inner)-int(sites[i] in inner)
        directed[sites[i],sites[j]] = d,a,b
        directed[sites[j],sites[i]] = d,-a,-b
        if d > 372:
            assert diff == a+b
        if d > 724:
            assert diff == b
    p,a,b = (1,5,29),(13,23,13),(23,13,13)
    cycles = [
        ((p,a,(23,7,17),p), (724,372,632), (1,0)),
        ((p,b,(7,23,17),p), (804,372,504), (0,1)),
        ((p,b,(1,29,5),(13,13,23),(29,1,5),a,p), (804,804,724,724,804,724), (-1,1)),
    ]
    for vertices,lengths,windings in cycles:
        traversed = [directed[i,j] for i,j in zip(vertices[:-1],vertices[1:],strict=True)]
        chord,clock724,clock804 = zip(*traversed,strict=True)
        assert chord == lengths and (sum(clock724),sum(clock804)) == windings
        print("闭路平方弦长与两钟净绕数",lengths,windings)


def check_support_counts(shell,angles):
    edges = shell["edges"]
    chord = {(i,j):d for i,j,d,_ in edges}
    clock_edges = [edge for edge,d in chord.items() if d in (724,804)]
    third = (6,9,9,7,9,6,10,12,13,14,14,13)
    occupied = set()
    for (i,j),k in zip(clock_edges,third,strict=True):
        legs = [tuple(sorted((i,k))),tuple(sorted((j,k)))]
        assert all(chord[e] not in (724,804) for e in legs)
        triangle = {(i,j),*legs}
        assert len(triangle) == 3 and occupied.isdisjoint(triangle)
        occupied.update(triangle)
    assert len(clock_edges) == 12

    inner = set(shell["inner"])
    charge = [int(site in inner) for site in shell["sites"]]
    zero = (0,0)
    for cutoff in (372,724):
        # 规范固定后，所有长边的比为一。长边同一连通分量上的其余规范角必须相等。
        group = list(range(len(charge)))
        for i,j,d,_ in edges:
            if d > cutoff:
                a,b = group[i],group[j]
                group = [min(a,b) if v in (a,b) else v for v in group]
        counts = defaultdict(Counter)
        for i,j,d,_ in edges:
            a,b = int(d==724),int(d==804)+charge[i]-charge[j]
            left,right = group[i],group[j]
            if cutoff == 372:
                a,b = 0,a+b
            if left > right:
                left,right,a,b = right,left,-a,-b
            counts[left,right][a,b] += 1
        unchanged = sum(c[zero] for c in counts.values())
        if cutoff == 372:
            assert len(set(group)) == 1
            assert set(counts[0,0]) == {zero,(0,1),(0,-1)}
            assert len(edges)-unchanged == 36
            continue
        labels = set().union(*(set(c) for c in counts.values()))
        assert labels == {zero,(0,1),(0,-1),(1,-1),(-1,1)}
        # 确认有限角标签在真实相位区间上模 2π 两两不同。
        for (a,b),(c,d) in combinations(labels,2):
            phase = acb(0,(a-c)*angles[0]+(b-d)*angles[1]).exp()
            assert abs(phase.arg()) > arb(".09")
        relaxed = sum(c[zero] if i==j else max(c.values()) for (i,j),c in counts.items())
        gains = {(i,j):max(c.values())-c[zero] for (i,j),c in counts.items()
                 if i!=j and max(c.values())>c[zero]}
        witnesses = {(3,11):15,(2,12):13,(0,18):9}
        assert gains == {edge:1 for edge in witnesses}
        used = set(witnesses)
        for (i,j),k in witnesses.items():
            for edge in (tuple(sorted((i,k))),tuple(sorted((k,j)))):
                assert edge not in used
                used.add(edge)
                c = counts[edge]
                assert c[zero] > max((mult for label,mult in c.items() if label!=zero),default=0)
        assert len(edges)-relaxed == 39
        assert len(edges)-unchanged == len(edges)-relaxed+len(witnesses) == 42
    print("最少改边数：无半径约束的非零两钟、D=372 的非零共同移、D=724 的实际钟移",
          12,36,42)


def check_phases(shell, model):
    up_index,down_index = 2345666,964890
    up,down = [acb.zeta_zero(n).imag for n in (up_index,down_index)]
    z_up,z_down = shell["phases"](up),shell["phases"](down)
    for qu,qd,expected in [
        (1,1,("3.02898420016",".09686584449","2.93211835568")),
        (-4,2,("-1.24530969908","1.50356763700","-2.74887733608")),
    ]:
        shifted = []
        for chord in (724,804):
            ratio = z_up[chord]**qu/z_down[chord]**qd
            omega = shell["theta"](chord)-shell["theta"](244)
            from_height = acb(0,omega*(qu*up-qd*down)).exp()
            assert abs(ratio-from_height) < arb("1e-30")
            shifted.append(ratio)
        angles = [v.arg() for v in shifted]+[(shifted[0]/shifted[1]).arg()]
        assert all(abs(x-arb(y)) < arb("1e-10") for x,y in zip(angles,expected,strict=True))
        print("两部门所取荷数、两钟有向角差及角差之差",qu,qd,[v.str(13) for v in angles])
        if (qu,qd) == (1,1):
            check_support_counts(shell,angles[:2])
    # 夸克窗中制备盒的中心是分别乘上部门荷数后的有效钟角。
    for charge,phases,centers in [
        (-4,z_up,("-.680061861483031",".873223965353125")),
        (2,z_down,(".565247837594779","-.630343671648397")),
    ]:
        for chord,center in zip((724,804),centers,strict=True):
            deviation = phases[chord]**charge/acb(0,arb(center)).exp()
            assert abs(deviation.arg()) < model["radius"]


if __name__ == "__main__":
    shell = run_path(str(Path(__file__).with_name("ch39-arc-relative-determinant.py")))
    model = run_path(str(Path(__file__).with_name("ch39-two-zero-quark-window.py")))
    check_graph(shell)
    check_phases(shell,model)
