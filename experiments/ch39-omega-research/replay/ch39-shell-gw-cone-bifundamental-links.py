"""Separate three-dimensional left and right scalar transports suffice for the measured flavor pair."""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb
import numpy as np

common = run_path(str(Path(__file__).with_name("ch39-shell-gw-cone-common-higgs.py")))
rank,native = common["rank"],common["native"]
q = arb


def certify():
    shell,full = rank["shell"],rank["full"]
    single,_,phases,_ = shell["transported_spectra"]()
    left,_ = rank["fermion"](1)
    right,_ = rank["fermion"](-1)
    theta = native["relative_phase"](q(676))
    a,d = [],[]
    # Fix the unobservable global zero-mode phase before comparing concrete matrices.
    anchor = acb("-.16257088602547","-.03649426870275")
    for i in range(2):
        z = native["p"][i]*phases[i]*left[i,0]*right[i,0]
        if i == 0:
            calibration = anchor*z.conjugate()/abs(anchor*z)
        z *= calibration
        z = acb(z.real.mid(),z.imag.mid())
        a.append(z)
        v = z.conjugate()*acb(0,-(-1)**i*theta).exp()
        d.append(acb(v.real.mid(),v.imag.mid()))
    det = a[0]*d[1]-a[1]*d[0]
    V,Du,Dd = native["arc"]["m"]["m"]["orb"]["reference"]()
    up = Du*q(".8616")
    down = (V*Dd*q(".01401")).conjugate()
    target = ((d[1]*up-a[1]*down)/det,(-d[0]*up+a[0]*down)/det)
    radius = q("1e-7")
    spatial = rank["normalized_trial"](single)
    R = q(12)/q(6).sqrt()
    for j,X in enumerate(target):
        adj = shell["adj"]
        gram = adj(X)*X
        singular = sorted(np.linalg.svd(full["model"]["array"](X),compute_uv=False))
        widths = q("1e-8")
        bounds = []
        for i,v in enumerate(singular):
            lo,hi = q(float(v))-widths,q(float(v))+widths
            assert full["count"](full["rows"](gram),lo**2) == i
            assert full["count"](full["rows"](gram),hi**2) == i+1
            bounds.append(lo.union(hi))
        shift = radius*(abs(a[1-j])+abs(d[1-j]))/abs(det)
        vals = [v+q(0,shift.upper()) for v in bounds]
        assert min(vals)>q(".0001") and vals[1]>q(".044") and vals[1]<q(".046")
        assert vals[-1]>13 and all(v<u for v,u in zip(vals,vals[1:]))
        assert sum(vals)<q("13.39")<3*R
    # With the ground vector's phase aligned, each true site differs in Frobenius
    # norm from an equally scaled rank-one matrix of norm R by < 2e-10.
    for i in range(6):
        assert abs(12*abs(spatial[i,0])-R)+q("1.21e-10") < q("2e-10")
    assert 12*q(2).sqrt()*(q("1e-9")+q("1e-12")) + 2*q(6).sqrt()*q("2e-10") < radius
    print("Output ball radius",radius.str(3),"odd/even target sums singular spectra",
          [list(np.linalg.svd(full["model"]["array"](X),compute_uv=False)) for X in target])


if __name__ == "__main__":
    certify()
