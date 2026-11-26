"""Positive flavor-singlet and traceless scattering witnesses below all 45 ultraviolet modes."""
from pathlib import Path
from runpy import run_path

from flint import arb

ten = run_path(str(Path(__file__).with_name("ch39-shell-gw-cone-ten-edges.py")))
q = arb


def certify():
    a,d = ten["coefficients"]()
    det = abs(a[0]*d[1]-a[1]*d[0])
    U = sum(q(z)**2 for z in (".00000615",".00311",".8616")).sqrt()
    D = sum(q(z)**2 for z in (".0000135",".000268",".01401")).sqrt()
    eps = q("1e-7")
    for i in range(2):
        assert (abs(d[i])*(q(".8616")-eps)-abs(a[i])*(D+eps))/det > q("13.28")
    # Fit imposes orthogonal lower bounds on the scaled low and high components.
    alpha = q(".4005869")/q(".408126")
    small = (U-4*q.pi()*q("1.001e-9"))/q(".408126")
    large = (alpha*U-D-4*q.pi()*q("1e-12"))/q(".078084")
    assert (small**2+large**2).sqrt()>q("10.858")
    # Frobenius/operator adjoint sum coefficient at the center. At all cone
    # points use the polar identification of the old and actual low bands.
    eta = q(".4081258")
    assert 6*abs(a[0])/q(6).sqrt() > eta
    assert max(abs(a[i]-a[0]) for i in range(2)) < q("1e-16")
    energy = q("1e5")
    mH2 = 3*q("60410.25")**2/q("1e-18")
    leakage = (4*q.pi())**2*energy**2/(mH2-energy**2)
    assert leakage < q("2e-16")
    for g, minimum in ((q("10.858"),"1.171"),(q(12),"1.43")):
        assert 3*((eta*g-eps)**2-leakage)/(16*q.pi())>q(minimum)
    for g, minimum in ((4*q.pi(),".597"),(q(12),".751")):
        gain = eta*(q("13.28")**2/g-g/2)-eps
        assert 3*(gain**2-leakage)/(16*q.pi()) > q(minimum)
    print("100 TeV frozen tree color-singlet J=0 lower bounds: 1.171 fully local;",
          ".597 even after eliminating the family-singlet Yukawa; strength 12 gives 1.43 and .751")


if __name__ == "__main__":
    certify()
