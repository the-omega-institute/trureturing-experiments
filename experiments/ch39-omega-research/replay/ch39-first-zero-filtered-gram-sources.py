"""同一有理滤波器作用于首个正零点的矩阵，制备五荷通用的三源处方。"""
from pathlib import Path
from runpy import run_path

from flint import acb_mat, arb, ctx


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    graph = lib("three-ray-chiral-boundary")
    ctx.prec = 256
    gram = lib("native-rh-ray-pvm")["gram_cholesky"](root)
    ports = [graph["soft"]["pole"],*graph["boundary_sites"]()]
    source = acb_mat([[gram[j][ports.index(i)] if i in ports else 0 for j in (1,2,3)]
                      for i in range(19)])
    assert sum(abs(z)**2 for z in source.entries()) < 4
    least = lib("terminal-port-light-bands")["least_singular"]
    spectrum = lib("adjacent-pvm-transport")["spectrum"]
    filter_error = max((arb("8.685")/9)**512, (9/arb("9.6"))**512)
    assert filter_error < arb("1.197e-8")
    for q in (1,-4,2,-3,6):
        levels,basis = spectrum(graph,q,1)
        assert 0 < levels[0] and levels[2] < arb("8.685") < 9 < arb("9.6") < levels[3]
        light = acb_mat([[basis[i,j] for j in range(3)] for i in range(19)])
        sigma = least(light.conjugate().transpose()*source)
        guaranteed = 7*sigma/9-2*filter_error
        assert sigma > arb(".0000176135") and guaranteed > arb(".0000136754")
        print("荷、第一枚正零点三源投影奇异值与任意目的背景滤波后下界",
              q,sigma.str(12),guaranteed.str(12))


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
