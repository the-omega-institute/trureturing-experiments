"""用低编号与两条弧的高零点极限约束单一纯端口的振幅与相位。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx

port = run_path(str(Path(__file__).with_name("ch39-single-port-ground-selection.py")))
study = port["study"]
native = study["native"]


def coordinates(T):
    return [T.trace().real/2, (T[0,0]-T[1,1]).real/2, T[0,1].real, T[0,1].imag]


def shifted_response(phase1,phase2):
    # Every candidate is evaluated at the same z in [0,2^-20].
    shell,block = study["phases_to_block"](phase1,phase2)
    diagonal = acb_mat(2,2)
    for i,k in enumerate(native["portals"]):
        diagonal[i,i] = shell[k,native["origin"]]
    J = native["source"]*diagonal*native["columns"].transpose()
    X = block.solve(J)
    adjoint = port["adjoint"]
    return coordinates(-adjoint(J)*X),coordinates(adjoint(X)*X)


def inequalities():
    height = acb.zeta_zero(11431).imag
    phase = study["geo"]["phases"](height)
    ref,derivative = shifted_response(phase[724],phase[804])
    zmax = arb(2)**-20
    shift = arb(0).union(zmax)
    gap = arb("1.3")-4*arb.pi()/32
    tail = 4*zmax**2/gap**3
    def difference(phase1,phase2):
        q,t = shifted_response(phase1,phase2)
        return [x-y+shift*(u-v)+arb(0,tail.upper()) for x,y,u,v in zip(q,ref,t,derivative)]

    low = []
    for n in (1,2,10,50,147):
        h = acb.zeta_zero(n).imag
        p = (h*h+study["offset"])/study["center"]
        assert 0 < p < 1
        phases = study["geo"]["phases"](h)
        low.append((p,difference(phases[724],phases[804])))
    rows = []
    for i,j in ((1,0),(-1,0),(0,1),(0,-1),(1,1),(-1,-1),(1,-1),(-1,1)):
        for k in (64,256,1024,4096):
            values = difference(phase[724]*acb(0,arb(i)/k).exp(),
                                phase[804]*acb(0,arb(j)/k).exp())
            for p,base in low:
                rows.append([p*l+(1-p)*r for l,r in zip(base,values)])
    return rows


def verify():
    ctx.prec = 192
    rows = inequalities()
    stack = [(0,0,1)]
    boxes = 0
    while stack:
        i,j,d = stack.pop()
        alpha = (arb(i)/d).union(arb(i+1)/d)*arb.pi()
        phi = (arb(j)/d).union(arb(j+1)/d)*arb.pi()
        # Any admissible port has the positive relative phase from the
        # five-zero real-port exclusion. Both endpoints are included.
        if (alpha > arb("2.220") and alpha < arb("2.242") and
                phi > arb(".892") and phi < arb(".920")):
            continue
        sine = alpha.sin()
        x,y,z = alpha.cos(),sine*phi.cos(),-sine*phi.sin()
        boxes += 1
        if any(a+b*x+c*y+e*z < 0 for a,b,c,e in rows):
            continue
        assert d < 2**19, (i,j,d)
        stack.extend((2*i+m,2*j+n,2*d) for m in (0,1) for n in (0,1))
    ratio = (arb("2.220").union(arb("2.242"))/2).tan()
    assert ratio > arb("2.014") and ratio < arb("2.072")
    print("Certified angular boxes",boxes,"amplitude ratio",ratio.str(5))


if __name__=="__main__":
    verify()
