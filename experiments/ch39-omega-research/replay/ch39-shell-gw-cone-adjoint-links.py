"""Measured flavor from a single SU(3) conjugation per shell site, with a verified submersion."""
import json
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx

folder = Path(__file__).resolve().parent
common = run_path(str(folder/"ch39-shell-gw-cone-common-higgs.py"))
rank,native = common["rank"],common["native"]
shell,full,q = rank["shell"],rank["full"],arb
adj,hs = shell["adj"],shell["hs"]


def seed_orbit(data,scale):
    cols = [acb_mat([[acb(*row[k])] for row in data]) for k in (0,1)]
    x = cols[0]/hs(cols[0],cols[0]).real.sqrt()
    y = cols[1]-x*hs(x,cols[1])
    y = y/hs(y,y).real.sqrt()
    return x*adj(y)*scale


def setup():
    _,_,phases,basis = shell["transported_spectra"]()
    l,_ = rank["fermion"](1)
    r,_ = rank["fermion"](-1)
    theta = native["relative_phase"](q(676))
    a,d = [],[]
    # Fix the unobservable global zero-mode phase before comparing concrete matrices.
    anchor = acb("-.16257088602547","-.03649426870275")
    for i in range(2):
        v = native["p"][i]*phases[i]*l[i,0]*r[i,0]
        if i == 0:
            calibration = anchor*v.conjugate()/abs(anchor*v)
        v *= calibration
        a.append(acb(v.real.mid(),v.imag.mid()))
        v = v.conjugate()*acb(0,-(-1)**i*theta).exp()
        d.append(acb(v.real.mid(),v.imag.mid()))
    V,Du,Dd = native["arc"]["m"]["m"]["orb"]["reference"]()
    P = acb_mat([[int((j-i)%3==1) for j in range(3)] for i in range(3)])
    targets = ((Du*P).transpose()*q(".8616"),(V*Dd*P*adj(V)).conjugate()*q(".01401"))
    seeds = json.loads((folder/"ch39-cone-conjugate-frames.json").read_text())
    return basis,a,d,targets,[[seed_orbit(v,12/q(6).sqrt()) for v in group] for group in seeds]


def certify():
    ctx.prec = 448
    basis,a,d,targets,sites = setup()
    zero = acb_mat(3,3)
    matrices = [sum((row[i]*sum(v,zero) for i,v in enumerate(sites)),zero) for row in (a,d)]
    assert sum(hs(x-y,x-y).real for x,y in zip(matrices,targets)).sqrt() < q("1e-10")
    U,D = matrices
    fields = [acb(0,1)*v for v in basis]
    basis = [shell["yuk"]["identity"](3)/q(3).sqrt()]+basis
    def real_coords(X):
        coeff = [hs(h,X) for h in basis]
        return [v.real for v in coeff]+[v.imag for v in coeff]
    cols = []
    for i,group in enumerate(sites):
        for H in group:
            for B in fields:
                variation = B*H-H*B
                cols.append(real_coords(a[i]*variation)+real_coords(d[i]*variation))
    for B in fields:
        cols.append(real_coords(B*U)+[q(0)]*18)
    for B in fields:
        cols.append([q(0)]*18+real_coords(D*B))
    J = acb_mat(cols).transpose()
    assert full["count"](full["rows"](J*adj(J)),q(1)/180**2) == 0
    # Unit-norm parameter direction. Each site conjugation plus one output
    # rotation has sum of generator norms <= sqrt(5) by Cauchy-Schwarz.
    bound = 5*(12/q(6).sqrt())*(sum((3*sum(abs(z) for z in row))**2 for row in (a,d))).sqrt()
    assert bound < 36
    radius = q("1e-7")
    model_error = q("1e-10")+12*q(2).sqrt()*(q("1e-9")+q("1e-12"))+2*q(6).sqrt()*q("2e-10")
    assert model_error+18*(180*radius)**2+q("5e-8") < radius
    print("36 real output coordinates; right inverse norm below 180; central physical target and",
          "output perturbations of Frobenius norm <=5e-8 achieved by SU(3) site conjugations")


if __name__ == "__main__":
    certify()
