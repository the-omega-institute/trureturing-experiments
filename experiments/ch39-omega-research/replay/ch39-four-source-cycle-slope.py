"""固定十项夸克读数的相位轨道上，四个二次源的壳三循环读数均非驻定。"""
from pathlib import Path
from runpy import run_path
import sys

from flint import acb, arb, arb_mat


def verify(replay):
    run_path(str(replay/"ch39-four-source-rotation-probe.py"))["verify"]()
    neighbours = run_path(str(replay/"ch39-four-quadratic-source-neighbours.py"))
    load = neighbours["load"]
    orbit, four = load("two-zero-orbit-condition"), load("two-zero-quadratic-source")
    model, _, _ = orbit["calibration"]()
    adj, trace, source = model["adj"], four["trace"], model["sources"]
    at_zero = load("quadratic-gram-source")["at_zero"]
    up, down = [at_zero(model,q,n) for q,n in ((-4,2345666),(2,964890))]
    parts, _ = orbit["orbit"](model,up,down)
    grams = [adj(m)*m for m in (up,down)]
    frame, length = four["normal_frame"](*(g/trace(g).real for g in grams))
    coefficients = four["orbit_coefficients"](parts,adj,frame)
    centers, radius = four["root_boxes"](coefficients,length)
    equations = arb_mat([[coefficients[i,j] for j in range(7)] for i in (1,2)])

    sites = model["shell"]["geo"]["sites"]
    cycle = [sites.index((y,z,x)) for x,y,z in sites]
    assert sorted(cycle) == list(range(len(sites)))
    assert all(cycle[cycle[cycle[i]]] == i for i in range(len(sites)))
    values = ("3.077473","3.273211","3.240416","3.161874")
    bounds = ((".10564",".10566"),("-.12241","-.12238"),("-.04360","-.04357"),(".06687",".06690"))
    for center, expected, (lo,hi) in zip(centers,values,bounds,strict=True):
        x,y = neighbours["refine"](center,radius,equations,four["trigonometric_columns"])
        phases = orbit["phases_at"](x,y)
        m = source*orbit["phase_sum"](parts,phases)
        dy = source*(acb(0,1)*phases[1]*parts[1])
        rho = m*adj(m)
        tangent = dy*adj(m) + m*adj(dy)
        # R e_i=e_{r(i)}，Tr(R rho)=sum_i rho[i,r(i)]。
        value = -2*sum(rho[i,j].imag for i,j in enumerate(cycle))
        derivative = -2*sum(tangent[i,j].imag for i,j in enumerate(cycle))
        assert abs(value-arb(expected)) < arb(".00001")
        assert arb(lo) < derivative < arb(hi)
        print("壳三循环读数的一阶相位导数",derivative.str(14))


if __name__ == "__main__":
    verify(Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent)
