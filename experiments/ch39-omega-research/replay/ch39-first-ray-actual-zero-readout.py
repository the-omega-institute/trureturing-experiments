"""直接在十个真实正零点重算三条射线首站对最低三态的最小奇异值。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, ctx


def verify(root):
    graph = run_path(str(root / "ch39-three-ray-chiral-boundary.py"))
    spectrum = run_path(str(root / "ch39-adjacent-pvm-transport.py"))["spectrum"]
    sites = graph["boundary_sites"]()
    assert sorted(sites) == [6, 7, 8]
    cases = ((1, 3992, 5351), (-4, 2958, 7359), (2, 1281, 3909),
             (-3, 3049, 4572), (6, 3837, 4623),
             (1, 22439, 807), (-4, 36298, 696), (2, 502189, 704),
             (-3, 19119, 629), (6, 60165, 744))
    for charge, number, millionths in cases:
        levels, basis = spectrum(graph, charge, number)
        assert levels[2] < 9 < levels[3]
        rows = acb_mat([[basis[i,j] for j in range(3)] for i in sites])
        values = (rows.conjugate().transpose()*rows).eig(algorithm="rump")
        assert all(v.is_finite() and v.imag.contains(0) and v.real > 0 for v in values)
        singular = min(v.real.sqrt().upper() for v in values)
        lower = min(v.real.sqrt().lower() for v in values)
        assert arb(millionths)/10**6 < lower <= singular < arb(millionths+1)/10**6
        print("荷、正零点编号、三射线首站最小奇异值", charge, number, singular.str(15))
    # h=ε/5，软辐边对算子漂移的贡献最终小于 ε/16。在 0<ε≤1 上 e/(g-e)<ε。
    ctx.prec = 144
    assert 2*(2*arb(2).sqrt()/5 + arb(1)/16) < arb("1.38785")
    assert (arb(".01")/(5*arb.pi()))**2 > arb("4.05e-7")


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
