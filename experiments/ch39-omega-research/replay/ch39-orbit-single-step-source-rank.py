"""两条最远旋转轨道的同背景单步秩一反例与固定首壳响应界。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


ROOT = Path(__file__).resolve().parent


def verify():
    def load(name):
        return run_path(str(ROOT / ("ch39-" + name + ".py")))

    graph = load("three-ray-chiral-boundary")
    pump = load("native-light-projector-pump")
    ctx.prec = 192
    ports = [14, 8, 7, 6]
    source_rows = load("native-rh-ray-pvm")["gram_cholesky"](ROOT)[1:]
    assert abs(source_rows[0][0]) > arb(".99")
    assert sum(abs(z) for z in source_rows[0][1:]) < arb(".041")
    source = acb_mat([[row[ports.index(i)] if i in ports else 0
                      for row in source_rows] for i in range(19)])
    least = load("terminal-port-light-bands")["least_singular"]
    edge = {(i,j):d for i,j,d,_ in graph["geo"]["edges"]}
    height = acb.zeta_zero(1).imag
    worst = arb(1)
    for charge in (1,-4,2,-3,6):
        base, (a724,a804) = pump["clock_matrices"](graph,charge)
        for sites, pivot in (((0,3,4),0),((1,2,5),1)):
            for i in sites:
                assert sorted(edge[i,j] for j in ports[1:]) == [244,724,804]
            assert all(base[i,14].contains(-1) for i in sites)
            coefficients = {edge[pivot,j]: (base+a724+a804)[pivot,j] for j in ports[1:]}
            a = coefficients[244]
            u,v = a / coefficients[724], a / coefficients[804]
            aligned = base+u*a724+v*a804+pump["adjoint"](u*a724+v*a804)
            assert all(aligned[i,j].contains(a) for i in sites for j in ports[1:])
            assert (abs(u)-1).contains(0) and (abs(v)-1).contains(0)
            if charge == 1:
                print("orbit, original phase angles", sites,u.arg().str(14),v.arg().str(14))
        rows = graph["laplacian"](charge,height)*source
        d = least(acb_mat([rows.tolist()[i] for i in (0,3,4)]))
        assert d > arb(".00050798")
        worst=min(worst,d)
        print("q, first-zero strict bound",charge,d.str(16))
    m=arb(".31895")
    d=arb(".00050798")
    g=arb("1.38785")
    s=(arb(".03")/g).sqrt()
    F=m*(1-s*s).sqrt()-s*(1-m*m).sqrt()
    assert m*m*d>arb(".00005167") and F*F*d>arb(".0000157")
    print("d*m^2, F, d*F^2", (d*m*m).str(15),F.str(15),(d*F*F).str(15))


if __name__ == "__main__":
    verify()
