"""以整数系数 (2,3,-64) 检验纯 11431 与规范均值的联合辨识。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx

resolution = run_path(str(Path(__file__).with_name("ch39-exposed-response-resolution.py")))
torus = resolution["certificate"]
A, B = arb(2)/64, arb(3)/64
curvature = arb(1)/5000
allowance = arb("7.5e-9")


def score(t):
    return A*(t[0,0]-t[1,1]).real/2+B*t[0,1].real-t[0,1].imag


def verify():
    ctx.prec = 192
    nodes,rows = torus["response_grid"]()
    det0,center = rows[0][0]
    reference = score(center)
    samples = []
    for i,row in enumerate(rows):
        values = []
        for j,(det,T) in enumerate(row):
            chord = 4-2*nodes[i].real-2*nodes[j].real
            values.append(det/det0*(score(T)-reference-curvature*chord+allowance))
        samples.append(values)
    torus["certify"](torus["torus_coefficients"](samples), stationary_origin=False)
    study = torus["study"]
    height = acb.zeta_zero(11431).imag
    h2 = acb.zeta_zero(2).imag
    c = study["offset"]
    C = height**2+c
    # The second zero crosses the center inside [1,2]; treat its affine
    # energy difference exactly, and exclude all the other competitors at once.
    resolution["spectrum_gap"](A,B,2/C,curvature,allowance,individual=(2,))
    other_gap = 1/C-1/(arb(resolution["cutoff"])**2+c)-allowance
    diff2 = score(study["native"]["plane_response"](h2)) - reference
    slope2 = 1/(h2*h2+c)-1/C
    crossing = diff2/slope2
    assert arb("1.97066597897") < crossing < arb("1.97066597899")
    assert diff2-slope2 > other_gap > arb("6.30e-10")
    print("Uniform off-center, non-2 gap on [1,2]; second-zero crossing",
          other_gap.str(14),crossing.str(18))
    # A concrete torus point bounds the strength required at the low end.
    u,v = arb(-161)/2**20, arb(544)/2**20
    phases = study["geo"]["phases"](height)
    matrix,block = study["phases_to_block"](phases[724]*acb(0,u).exp(),phases[804]*acb(0,v).exp())
    native = study["native"]
    solution = block.solve(native["source"])
    ports = native["portals"]
    origin = native["origin"]
    small = acb_mat([[-matrix[origin,i]*solution[native["inside"].index(i),k]*matrix[j,origin]
                      for k,j in enumerate(ports)] for i in ports])
    probe = native["columns"]*small*native["columns"].transpose()
    obstruction = -C*(score(probe)-reference)
    assert obstruction > arb(".8838")
    print("High-zero obstruction for all strengths up to",obstruction.str(15))


if __name__ == "__main__":
    verify()
