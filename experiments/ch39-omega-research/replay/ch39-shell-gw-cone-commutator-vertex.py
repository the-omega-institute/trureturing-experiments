"""A local double commutator and the two quark scattering blocks above the eight light thresholds."""
from pathlib import Path
from runpy import run_path

from flint import arb

folder = Path(__file__).resolve().parent
ten = run_path(str(folder/"ch39-shell-gw-cone-ten-edges.py"))
q = arb


def certify():
    a,d = ten["coefficients"]()
    for row in (a,d):
        assert sum(abs(z)**2 for z in row)*3 < q(".408126")**2
    theta = ten["orbit"]["native"]["relative_phase"](q(676))
    up,down = [sum(q(v)**2 for v in group).sqrt() for group in (
        (".00000615",".00311",".8616"),(".0000135",".000268",".01401"))]
    error = q("1e-7")
    # Ground-normalized D² fixes both extreme root matrices and acts by 1/4
    # on the middle root spaces. Square again when forming the exchange Gram.
    full_norm = 12*q(".408126")+error
    extremal = [2*(up+error)**2,
                (down+error)**2+(2*abs(theta.cos())*up+down+error)**2]
    positive = []
    light_weight = 1/(1-(q("25618")/q("1e5"))**2)
    for c,limit in zip(extremal,(".185",".271")):
        light_gram = full_norm**2/16+15*c/16+error
        bound = 3*light_weight*light_gram/(16*q.pi())
        assert bound < q(limit)
        positive.append(bound)
    assert sum(positive) < q(".455")
    Q, _ = ten["rank"]["cone"]["first"]["gauge_target"]()
    def negative(t):
        heavy_mass2 = 3*q("60410.25")**2/t
        return 3*full_norm**2*Q**2/((heavy_mass2-Q**2)*16*q.pi())
    assert negative(q("1e-18")) < q(".425")
    assert 2*negative(q("6e-19")) < q(".46")
    # A contractive filter changes the light-ground readout by at most twice
    # the global deviation from its limiting extreme eigenvector.
    assert 2e-13+12*q(2).sqrt()*(q("2.1e-15")+q("4e-15")+q("1e-17"))+15*(5200*q("5e-10"))**2+q("3e-10") < q("5e-10")
    print("positive quark block bounds",*[v.upper().str(8) for v in positive],
          "negative bound through Q",negative(q("1e-18")).upper().str(8))


if __name__ == "__main__":
    certify()
