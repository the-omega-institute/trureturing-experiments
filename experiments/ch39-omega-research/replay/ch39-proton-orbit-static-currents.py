"""三循环荷五静态传播的全域正下界及八站上的实际零点相消。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx


ROOT = Path(__file__).resolve().parent


def verify():
    graph = run_path(str(ROOT / "ch39-three-ray-chiral-boundary.py"))
    clock = run_path(str(ROOT / "ch39-native-light-projector-pump.py"))
    ctx.prec = 192
    sites = (0, 2, 3, 5, 6, 9, 11, 13)
    def green(n):
        return graph["laplacian"](5, acb.zeta_zero(n).imag).inv()

    propagator = green(177)
    x, y, z = [propagator[6,j] for j in (3, 5, 13)]
    areas = [(u.conjugate()*v).imag for u, v in ((y,z), (z,x), (x,y))]
    assert all(t < 0 for t in areas)
    print("第 177 枚零点三个非负电流权", [w.str(12) for w in (t/sum(areas) for t in areas)])
    first = green(1)
    assert all(first[i,j].real > arb(".006") for i in sites for j in sites)

    ctx.prec = 96
    base, links = clock["clock_matrices"](graph, 5)
    adjoint = clock["adjoint"]
    n = 64
    roots = [acb(0, 2*arb.pi()*k/n).exp() for k in range(n)]
    layers = [[link*u + adjoint(link*u) for u in roots] for link in links]
    error = 4*arb(2).sqrt()*(arb.pi()/(2*n)).sin() + 6*arb("1e-8")
    def squared_abs(u):
        return u.real**2 + u.imag**2

    for left in layers[0]:
        for right in layers[1]:
            center = (base+left+right).inv()
            denominator = 1-error*sum(map(squared_abs, center.entries())).sqrt()
            assert denominator > 0
            column_squared = sum(squared_abs(center[k,0]) for k in range(19))
            assert center[0,3].real-error*column_squared/denominator > arb(".0038154")
    print("三站轨道上的非负单位电流下界", ".0038154")


if __name__ == "__main__":
    verify()
