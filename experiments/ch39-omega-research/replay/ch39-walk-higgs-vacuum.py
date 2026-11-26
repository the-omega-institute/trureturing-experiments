"""荷三有限扩散步数的任意非负混合都保持真空顶底比的统一上界。"""
from pathlib import Path
from runpy import run_path
import json

from flint import acb,acb_mat,arb,ctx


def verify():
    root = Path(__file__).parent
    shell = run_path(str(root/"ch39-three-ray-chiral-boundary.py"))
    flavor = run_path(str(root/"ch39-minimum-flavor-support.py"))
    vacuum = run_path(str(root/"ch39-higgs-phase-vacuum.py"))
    positive = run_path(str(root/"ch39-three-mode-higgs-hierarchy.py"))["positive"]
    ctx.prec = 320
    adjoint = flavor["adjoint"]
    up,down = flavor["flavor_maps"](shell,11431)
    penalty = 1000*adjoint(down)*down-195*adjoint(up)*up
    lap = shell["laplacian"](3,acb.zeta_zero(11431).imag)
    values, modes = lap.eig(right=True)
    order = sorted(range(19),key=lambda j:float(values[j].real))
    levels = [values[j].real for j in order]
    assert all(a < b for a,b in zip(levels,levels[1:])) and levels[-1] < 36
    scale = 36-levels[0]
    sites = shell["geo"]["sites"]
    rot = [sites.index((y,z,x)) for x,y,z in sites]
    data = json.loads((root/"ch39-phase-vacuum-witnesses.json").read_text())
    w = vacuum["gram"](rot,data["hessian"])
    positive(w)
    distances = acb_mat([[w[i,i]+w[j,j]-2*w[i,j] for j in range(19)] for i in range(19)])
    def tested(matrix):
        return penalty + acb_mat([[distances[i,j]*matrix[i,j] for j in range(19)] for i in range(19)])
    identity = acb_mat([[int(i == j) for j in range(19)] for i in range(19)])
    step = (36*identity-lap)/scale
    current = identity
    for n in range(1,100):
        current = current*step
        positive(tested(-scale*current))
    ground = acb_mat([[modes[i,order[0]]] for i in range(19)])
    projector = (ground*adjoint(ground))/(adjoint(ground)*ground)[0,0]
    margin = arb(".001")
    positive(tested(-scale*projector)-margin*identity)
    ratio = (36-levels[1])/scale
    assert ratio < arb(".818")
    max_entry = max(abs(z) for row in distances.tolist() for z in row)
    assert 19*max_entry*scale*ratio**100 < margin
    print("核查的正步数、余下步数起点、次模衰减比、任意非负混合的真空顶底比上界",
          99,100,ratio.str(10),(arb(3000)/195).sqrt().str(12))


if __name__ == "__main__":
    verify()
