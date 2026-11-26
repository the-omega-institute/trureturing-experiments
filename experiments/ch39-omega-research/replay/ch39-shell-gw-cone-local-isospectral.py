"""Local unitary scalar transport can align the actual light line to a common fitted Higgs."""
from pathlib import Path
from runpy import run_path

from flint import acb, arb

common = run_path(str(Path(__file__).with_name("ch39-shell-gw-cone-common-higgs.py")))
rank,native = common["rank"],common["native"]
q = arb


def certify():
    shell = rank["shell"]
    single,_,phases,_ = shell["transported_spectra"]()
    spatial = rank["normalized_trial"](single)
    left,_ = rank["fermion"](1)
    right,_ = rank["fermion"](-1)
    phase = native["relative_phase"](q(676))
    reference = native["p"]
    a = [reference[i]*phases[i]*left[i,0]*right[i,0] for i in range(6)]
    d = [z.conjugate()*acb(0,-(-1)**i*phase).exp() for i,z in enumerate(a)]
    # Use a concrete, exactly specified row for the even/odd model and
    # certify its distance from the full true rows, including vector residuals.
    A = [acb(a[i].real.mid(),a[i].imag.mid()) for i in range(2)]
    D = [acb(d[i].real.mid(),d[i].imag.mid()) for i in range(2)]
    for actual,base in ((a,A),(d,D)):
        assert common["norm"]([z-base[i%2] for i,z in enumerate(actual)])+q("2e-11") < q("1e-9")
    radii = [12*abs(spatial[i,0]) for i in range(6)]
    assert all(q("4.89897")+q("2e-10") < z < q("4.89899")-q("2e-10") for z in radii)
    det = abs(A[0]*D[1]-A[1]*D[0])
    U = common["norm"]([q(x) for x in (".00000615",".00311",".8616")])
    V = common["norm"]([q(x) for x in (".0000135",".000268",".01401")])
    radius = q("1e-5")
    for i in range(2):
        # Bound the norm of both required three-site matrix sums.
        upper = (abs(D[i])*(U+radius)+abs(A[i])*(V+radius))/det
        lower = (abs(D[i])*(U-radius)-abs(A[i])*(V+radius))/det
        assert q("13.29") < lower and upper < q("13.75")
        variation = (abs(D[i])+abs(A[i]))*radius/det
        assert variation < q(".001")
    # Transverse radii in each three-site group permit strict triangle closure.
    min_perp = (q("4.89897")**2-(q("13.75")/3)**2).sqrt()
    max_perp = (q("4.89899")**2-(q("13.29")/3)**2).sqrt()
    assert max_perp < 2*min_perp
    assert 12*q(2).sqrt()*(q("1e-9")+q("1e-12")) < radius
    print("Strength 12, target ball radius",radius.str(3),"six allowed site lengths (4.89897,4.89899)",
          "individual transverse radii", min_perp.lower().str(7),max_perp.upper().str(7))


if __name__ == "__main__":
    certify()
