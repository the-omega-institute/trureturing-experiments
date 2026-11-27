"""由十九点壳的闭路与顶点重定相验证单位荷两钟联络的两个局域阈值。"""
from fractions import Fraction
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
