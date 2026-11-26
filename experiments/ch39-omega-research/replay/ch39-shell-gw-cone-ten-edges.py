"""Ten of the twelve shell edges suffice at common strength 12 throughout the eight-threshold cone."""
import json
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx

folder = Path(__file__).resolve().parent
orbit = run_path(str(folder/"ch39-shell-gw-cone-adjoint-links.py"))
rank,shell,full = orbit["rank"],orbit["shell"],orbit["full"]
hs,adj,eye = shell["hs"],shell["adj"],shell["yuk"]["identity"]
q = arb


def refine(K,shift,gap,zero=False):
    v = rank["normalized_trial"](K)
    for _ in range(2):
        v = (K-shift*eye(K.nrows())).solve(v)
        v /= hs(v,v).real.sqrt()
    e = q(0) if zero else hs(v,K*v).real
    rem = K*v-e*v
    assert 2*hs(rem,rem).real.sqrt()/gap < q("1e-19")
    return v


def refined_fermion(sign):
    K = shell["yuk"]["model"]["action_kernel"](sign,q(0),q(0))[1]
    K = acb_mat([[K[3*i][3*j] for j in range(12)] for i in range(12)])
    b,x = q(".33")/8,K/8
    unit = eye(12)
    while 1-b > q("1e-36"):
        x = x*(3*unit-x*x)/2
        b = b*(3-b*b)/2
    o = 0 if sign>0 else 6
    block = acb_mat([[(unit[i,j]-sign*x[i,j])/2 for j in range(o,o+6)] for i in range(o,o+6)])
    assert full["count"](full["rows"](block),q(".0051")) == 1
    return refine(block,q("1e-12"),q(".005"),True)


def coefficients():
    single,_,phases,_ = shell["transported_spectra"]()
    v = rank["normalized_trial"](single)
    e = hs(v,single*v).real
    v = refine(single,e+q("1e-12"),q("3.6"))
    ref = orbit["native"]["p"]
    ref0 = v[0,0]/ref[0]
    for i in range(6):
        assert abs(v[i,0]/ref[i]-ref0/abs(ref0)/q(6).sqrt()) < q("1e-18")
    l,r = refined_fermion(1),refined_fermion(-1)
    a = [ref[i]*phases[i]*l[i,0]*r[i,0] for i in range(6)]
    anchor = acb("-.16257088602547","-.03649426870275")
    calibration = anchor*a[0].conjugate()/abs(anchor*a[0])
    a = [z*calibration for z in a]
    theta = orbit["native"]["relative_phase"](q(676))
    d = [z.conjugate()*acb(0,-(-1)**i*theta).exp() for i,z in enumerate(a)]
    for row in (a,d):
        assert all(abs(z-row[i%2]) < q("1e-18") for i,z in enumerate(row))
    return a[:2],d[:2]


def frame(data):
    Q = []
    for k in range(3):
        v = acb_mat([[acb(*row[k])] for row in data])
        for w in Q:
            v -= w*hs(w,v)
        Q.append(v/hs(v,v).real.sqrt())
    return acb_mat([[v[i,0] for v in Q] for i in range(3)])


def certify():
    ctx.prec = 480
    a,d = coefficients()
    g,_,_,targets,_ = orbit["setup"]()
    raw = json.loads((folder/"ch39-cone-ten-edge-frames.json").read_text())
    Fs = [frame(f) for f in raw]
    w,u,y = Fs[-3:]
    mats = [orbit["seed_orbit"](raw[i],(2 if i%2==0 else 1)*12/q(6).sqrt()) for i in range(4)]
    h = [sum((row[i//2]*m for i,m in enumerate(mats)),acb_mat(3,3)) for row in (a,d)]
    up = u*targets[0]*w.transpose()
    down = w.conjugate()*targets[1]*y
    assert (hs(up-h[0],up-h[0])+hs(down-h[1],down-h[1])).real.sqrt() < q("2e-13")
    basis = [eye(3)/q(3).sqrt()]+g
    B = [acb(0,1)*b for b in basis]
    def coords(X):
        coeff = [hs(b,X) for b in basis]
        return [z.real for z in coeff]+[z.imag for z in coeff]
    cols = []
    for i,H in enumerate(mats):
        for b in B[1:]:
            z = b*H-H*b
            cols.append(coords(a[i//2]*z)+coords(d[i//2]*z))
    for b in B:
        cols.append(coords(-up*b.transpose())+coords(-b.conjugate()*down))
    for b in B:
        cols.append(coords(-b*up)+[0]*18)
    for b in B:
        cols.append([0]*18+coords(-w.conjugate()*targets[1]*b*y))
    J = acb_mat(cols).transpose()
    assert full["count"](full["rows"](J*adj(J)),1/q(5200)**2)==0
    R = 12/q(6).sqrt()
    shell_second = 4*R*sum((3*sum(abs(z) for z in row))**2 for row in (a,d)).sqrt()
    physical_second = 2*sum(hs(X,X).real for X in targets).sqrt()
    assert shell_second+physical_second < 30
    # Tight residuals above and the previously certified cone field/fermion motion.
    # The unit-norm scalar site distortion here is < 7e-16 globally,
    # before the 12-fold rescaling; the two fermion frame motions sum to <4e-15.
    t = q("1e-18")
    D = 8+q(16)/3
    assert 2*(4*t/3+D*t/(q(".08991275")/2-D*t)) < q("7e-16")
    error = q("2e-13")+12*q(2).sqrt()*(q("7e-16")+q("4e-15")+q("1e-17"))
    radius = q("5e-10")
    assert error+15*(5200*radius)**2+q("3e-10") < radius
    print("four shared site conjugations suffice; two exact unaltered long edges, ten necessary and sufficient",
          "right inverse norm below 5200, perturbation in target norm up to 3e-10")


if __name__ == "__main__":
    certify()
