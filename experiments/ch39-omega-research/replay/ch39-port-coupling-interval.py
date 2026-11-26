"""七三角形单端口的连续耦合段与首次全谱失选的范围。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx

point = run_path(str(Path(__file__).with_name("ch39-port-coupling-competition.py")))
torus, resolution = point["torus"],point["resolution"]
study, native = point["study"],point["native"]
adjoint = point["adjoint"]
size, degree = 80, 38


def block_and_source(first,second,vector):
    shell,D = study["phases_to_block"](first,second)
    m = acb_mat(2,2)
    for i,k in enumerate(native["portals"]):
        m[i,i] = shell[k,native["origin"]]
    b = native["source"]*m*native["columns"].transpose()*vector
    return D,b


def controls(D,b,lo,hi):
    data = []
    for z,sign in [(lo,1),(hi,-1)]:
        shifted = acb_mat(D)
        for i in range(D.nrows()):
            shifted[i,i] += z
        x = shifted.solve(b)
        f = -(adjoint(b)*x)[0,0].real
        deriv = sum(abs(t)**2 for t in x.entries())
        det = shifted.det().real
        assert det>0
        data.extend(((det**2,f),(det**2,f+sign*(hi-lo)*deriv/3)))
    return data


def fourier(grid):
    first = [acb.dft(row) for row in grid]
    cols = [acb.dft(col) for col in zip(*first)]
    return {(i,j):cols[j%size][i%size]/size**2 for i in range(-degree,degree+1)
            for j in range(-degree,degree+1)}


def certify_interval(lo,hi,budget,curvature,gain):
    vector = point["prepare"]()
    ref_h = acb.zeta_zero(11431).imag
    phase = study["geo"]["phases"](ref_h)
    ref = controls(*block_and_source(phase[724],phase[804],vector),lo,hi)
    nodes = [acb(0,2*arb.pi()*k/size).exp() for k in range(size)]
    grids = [[] for _ in ref]
    for i,u in enumerate(nodes):
        lines = [[] for _ in ref]
        for j,v in enumerate(nodes):
            vals = controls(*block_and_source(phase[724]*u,phase[804]*v,vector),lo,hi)
            chord = 4-2*nodes[i].real-2*nodes[j].real
            for line,(d,y),(d0,y0) in zip(lines,vals,ref):
                line.append(d/d0*(y-y0+budget-curvature*chord))
        for grid,line in zip(grids,lines):
            grid.append(line)
    for k,grid in enumerate(grids):
        print("Hermite control",k,flush=True)
        torus["certify"](fourier(grid),stationary_origin=False)
    gap_D = arb("1.3")-4*arb.pi()/32
    remainder = (hi-lo)**4/(4*(gap_D+lo)**5)
    def delta(h):
        if h.overlaps(ref_h):
            return arb(0)
        D,J = point["internal"](h)
        vals = [y-y0 for (_,y),(_,y0) in zip(controls(D,J*vector,lo,hi),ref)]
        low = min(y.lower() for y in vals)
        high = max(y.upper() for y in vals)
        return low.union(high)+arb(0,remainder.upper())
    gap,_ = resolution["spectrum_gap"](0,0,1/study["center"],curvature,budget+remainder,
                                       readout=lambda y:y,response_at=delta)
    z = point["negative_energy"](ref_h,vector,gain)
    assert z < hi
    # z(g) increases with g. Its g^2 ratio is at least the endpoint value.
    norm = 1+2*gain**2/gap_D**2
    multiple = min((z/gain**2).lower(),(gap/norm).lower())
    print("Hermite remainder, gap per squared coupling up to",gain.str(7),
          remainder.str(11),multiple.str(12))
    return multiple


if __name__=="__main__":
    ctx.prec = 192
    first = certify_interval(arb(0),arb(1)/200,arb("5e-9"),arb(1)/5000,arb(1)/4)
    assert first > arb("2.45e-9")
    left = arb(".00465")
    reference = acb.zeta_zero(11431).imag
    vector = point["prepare"]()
    assert point["negative_energy"](reference,vector,arb(1)/4) > left
    second = certify_interval(left,arb("0.007335"),arb("8.04e-9"),arb(1)/20000,arb(5)/16)
    assert second > arb("5e-11")

    competitor = acb.zeta_zero(1314155354).imag
    gain = arb(63)/200
    difference = (point["negative_energy"](competitor,vector,gain)
                  - point["negative_energy"](reference,vector,gain))
    assert difference > arb("5.87e-12")
    print("Coupling",gain.str(8),"zero 1314155354 ordinate and energy advantage",
          competitor.str(20),difference.str(15))
