"""同一 26 条路径七三角端口在中等与单位耦合下的全谱选点。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx

small = run_path(str(Path(__file__).with_name("ch39-single-port-ground-selection.py")))
study = small["study"]
native, geo = study["native"],study["geo"]
torus, resolution = small["torus"],small["resolution"]
adjoint = small["adjoint"]


def prepare():
    loops = [((5,13,14),1),((5,13,15),4),((5,13,17),1),((5,15,16),1)]
    phase = acb(1)
    for (i,j,k),num in loops:
        z = geo["unit"](i,j)*geo["unit"](j,k)*geo["unit"](i,k).conjugate()
        if z.imag < 0:
            z = z.conjugate()
        assert z.imag > 0
        phase *= z**num
    return acb_mat([[(arb(5)/26).sqrt()],[phase*(arb(21)/26).sqrt()]])


def internal(h):
    phase = geo["phases"](h)
    shell,D = study["phases_to_block"](phase[724],phase[804])
    diagonal = acb_mat(2,2)
    for i,j in enumerate(native["portals"]):
        diagonal[i,i] = shell[j,native["origin"]]
    return D,native["source"]*diagonal*native["columns"].transpose()


def solve(D,J,z):
    shifted = acb_mat(D)
    for i in range(D.nrows()):
        shifted[i,i] += z
    return shifted.solve(J)


def negative_energy(h,vector,gain, *, gauge_weight=1):
    D,J = internal(h)
    b = J*vector
    inv = gauge_weight/(h*h+study["offset"])
    def equation(z):
        return z-gain**2*(inv+(adjoint(b)*solve(D,b,z))[0,0].real)
    lo,hi = arb(0),arb(1)
    assert equation(lo)<0 and equation(hi)>0
    for _ in range(100):
        mid = (lo+hi)/2
        value = equation(mid)
        if value>0:
            hi = mid
        else:
            assert value<0
            lo = mid
    return lo.union(hi)


def full_spectrum(h,vector,gain, *, curvature=None, allowance=None, gauge_weight=1):
    shift = negative_energy(h,vector,gain,gauge_weight=gauge_weight)
    nodes,rows = torus["response_grid"](shift)
    readout = lambda t: (adjoint(vector)*t*vector)[0,0].real
    det0,T = rows[0][0]
    ref = readout(T)
    if curvature is None:
        curvature = arb(1)/5000
    if allowance is None:
        allowance = arb("5e-9")
    values = [[d/det0*(readout(t)-ref+allowance-curvature*(4-2*nodes[i].real-2*nodes[j].real))
               for j,(d,t) in enumerate(row)] for i,row in enumerate(rows)]
    torus["certify"](torus["torus_coefficients"](values),stationary_origin=False)
    def response(h):
        D,J = internal(h)
        return -adjoint(J)*solve(D,J,shift)
    gap,_ = resolution["spectrum_gap"](0,0,gauge_weight/study["center"],curvature,allowance,
                                       readout=readout,response_at=response)
    # Scalar Schur derivative magnitude <= 1+g^2 ||Jv||^2 / gap(D)^2
    # on the path between the negative roots.
    slope = 1+gain**2*2/(arb("1.3")-4*arb.pi()/32)**2
    full = min(shift.lower(),(gain**2*gap/slope).lower())
    assert full > 0
    print("Coupling",gain.str(8),"ground energy",-shift,
          "effective and full gap lower",gap.str(12),full.str(12))
    return full


def verify():
    ctx.prec = 192
    ref = acb.zeta_zero(11431).imag
    vector = prepare()
    assert full_spectrum(ref,vector,arb(1)/4) > arb("1.69e-10")
    h = acb.zeta_zero(1080815).imag
    gap = negative_energy(h,vector,arb(1))-negative_energy(ref,vector,arb(1))
    assert gap > arb("1e-7")
    print("Unit coupling, zero 1080815 energy below 11431 by",gap.str(15))


if __name__ == "__main__":
    verify()
