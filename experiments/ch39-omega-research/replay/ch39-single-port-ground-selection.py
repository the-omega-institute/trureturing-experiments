"""用一个高斯有理端口与原十八维壳内部块选出纯 11431。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx

resolution = run_path(str(Path(__file__).with_name("ch39-exposed-response-resolution.py")))
torus = resolution["certificate"]
study = torus["study"]
v = acb_mat([[acb(26,16)/64],[acb(7,62)/64]])
curvature = arb(1)/5000
allowance = arb("3e-9")


def adjoint(matrix):
    return matrix.conjugate().transpose()


def verify(vector=v, *, allowance=allowance):
    def readout(T):
        return (adjoint(vector)*T*vector)[0,0].real
    ctx.prec = 192
    nodes,rows = torus["response_grid"]()
    det0,center = rows[0][0]
    reference = readout(center)
    values = []
    for i,row in enumerate(rows):
        values.append([det/det0*(readout(T)-reference+allowance
                                 -curvature*(4-2*nodes[i].real-2*nodes[j].real))
                       for j,(det,T) in enumerate(row)])
    torus["certify"](torus["torus_coefficients"](values), stationary_origin=False)
    h = acb.zeta_zero(11431).imag
    c = study["offset"]
    C = h*h+c
    gap,_ = resolution["spectrum_gap"](0,0,1/C,curvature,allowance,readout=readout)
    assert gap > 0
    # Show that the original internal block and a single column realize this score.
    native = study["native"]
    phase = study["geo"]["phases"](h)
    shell,D = study["phases_to_block"](phase[724],phase[804])
    phase_matrix = acb_mat(2,2)
    for i,portal in enumerate(native["portals"]):
        phase_matrix[i,i] = shell[portal,native["origin"]]
    source = native["source"]*phase_matrix*native["columns"].transpose()*vector
    solved = D.solve(source)
    assert abs((adjoint(source)*solved)[0,0]+reference) < arb("1e-32")
    eps = arb(2)**-10
    center_norm = sum(abs(z)**2 for z in solved.entries())
    correction = eps**2*(1/C-reference)*center_norm
    full_gap = eps**2*(gap-correction)
    assert full_gap > 0
    print("Single-port effective gap",gap.str(13),"coupling",eps.str(10),
          "center correction in effective units",correction.str(9),
          "full nineteen-dimensional ground gap",full_gap.str(15))
    # Exact central root, bracketed by Schur evaluation at two real energies.
    upper = eps**2*(reference-1/C+correction)
    lower = eps**2*(reference-1/C)
    eye = acb_mat([[int(i==j) for j in range(18)] for i in range(18)])
    def schur(E):
        return -eps**2/C-E-eps**2*(adjoint(source)*(D-E*eye).solve(source))[0,0].real
    # Both theoretical inequalities are strict at the center.
    assert schur(lower.lower())>0 and schur(upper.upper())<0
    print("Ground energy enclosure",lower.union(upper).str(15))
    return gap,full_gap,lower.union(upper)


if __name__ == "__main__":
    gap,full_gap,_ = verify()
    assert gap > arb("5.12e-9") and full_gap > arb("2.15e-15")
