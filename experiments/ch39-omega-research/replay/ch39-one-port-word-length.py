"""实际第七壳荷六磁图上，末端口向量为循环向量。"""
from pathlib import Path
from runpy import run_path

from flint import arb


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))
    graph=lib("three-ray-chiral-boundary")
    port=graph["boundary_sites"]()[-1]
    assert port==6
    _,basis=lib("adjacent-pvm-transport")["spectrum"](graph,6,2_306_888_100_092)
    born=[abs(basis[port,j])**2 for j in range(19)]
    assert all(w.is_finite() and w>arb(".00008498") for w in born)
    print("第七壳十九条本征线在末端口的最小 Born 权下界",min(w.lower() for w in born).str(13),flush=True)


if __name__=="__main__":
    verify(Path(__file__).resolve().parent)
