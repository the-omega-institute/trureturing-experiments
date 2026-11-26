"""Continuum gauge feasibility, true-zero pinning, and a global normalized Hardy penalty."""
from fractions import Fraction
from pathlib import Path
from runpy import run_path
from flint import acb, arb, arb_series, ctx

gauge = run_path(str(Path(__file__).with_name("ch39-floating-scale-winner.py")))
q = gauge["q"]


def modulus(u):
    return (3-2*u)*(2*(113+72*u)).sqrt()/13


def weight(g):
    return modulus(arb(2).sqrt()*(g*arb(2).log()).cos())/(g*g+arb(1)/4)


def hardy(g):
    return arb_series([g]).riemann_siegel_z()[0]


def verify():
    ctx.prec, ctx.cap = 144, 4
    root = acb.zeta_zero(10875).imag
    best = weight(root)
    variable = arb_series([root,1])
    z, w = variable.riemann_siegel_z(), weight(variable)
    assert z[0].contains(0) and z[1]>5 and w[1]<0
    critical = -w[1]/(best*z[1])
    assert arb(".00141357") < critical < arb(".00141359")
    assert w[2]/best + critical*z[2] < 0
    offset = q(gauge["fixed"]["OFFSET"])
    common = gauge["target"](10875,arb("26.859482438"))
    delta = arb(".003")
    band = root+arb(0,delta)
    assert abs(common*((root**2+offset)/(band**2+offset)-1)) < arb(".000002")
    assert weight(root-delta)-best > arb("6.87e-13")

    # Extend the published simultaneous gauge exclusion to all heights below 10603.
    low_zero = acb.zeta_zero(10854).imag
    time = q(gauge["fixed"]["L"]+Fraction(".04860"))
    low_target = gauge["target"](10854,time)
    assert low_target*(1-(low_zero**2+offset)/(arb(10603)**2+offset)) < arb(".00049")
    cap = (modulus(-arb(59)/108)/best-arb(1)/4).sqrt()
    assert cap < arb("10621.28")
    # Enclose the derivative across the central cell. Subdivide to control Z' cancellation.
    near = Fraction(1,1024)
    for j in range(-32,32):
        g = root + q(near*(2*j+1)/64) + arb(0,q(near/64))
        v = arb_series([g,1])
        assert v.riemann_siegel_z()[1] > 5
        assert abs(weight(v)[1])/best < arb(".015")
    outside_ratio = arb(".003641")
    assert arb(".015") < arb(5)*outside_ratio

    # For g below the central cell and above it, prove W/W* - 1 < |Z|/256.
    left = (root-q(near)).upper()
    right = (root+q(near)).lower()
    peak_left, peak_right = arb("10617.639").lower(), arb("10617.652").upper()
    stack = [(arb(10603),peak_left),(peak_right,left),(right,arb("10621.28"))]
    accepted = 0
    while stack:
        lo,hi = stack.pop()
        cell = lo.union(hi)
        excess = weight(cell)/best-1
        if excess<0 or excess < outside_ratio*abs(hardy(cell)):
            accepted += 1
            continue
        assert hi-lo > arb("1e-9")
        middle = (lo+hi)/2
        stack.extend([(lo,middle),(middle,hi)])
    def penalty_ratio(g):
        variable = arb_series([g,1])
        z = variable.riemann_siegel_z()
        assert not z[0].contains(0)
        return (weight(variable)/best-1)/(z if z[0]>0 else -z)

    # Certify strict concavity throughout the only remaining interval.
    pieces = [(peak_left,peak_right)]
    while pieces:
        lo,hi = pieces.pop()
        cell = lo.union(hi)
        if not hardy(cell).contains(0) and penalty_ratio(cell)[2] < 0:
            continue
        assert hi-lo > arb("1e-8")
        mid = (lo+hi)/2
        pieces.extend([(lo,mid),(mid,hi)])
    lo,hi = arb("10617.645332791202380"), arb("10617.645332791202395")
    assert peak_left < lo < hi < peak_right
    assert penalty_ratio(lo)[1]>0>penalty_ratio(hi)[1]
    enclosure = lo.union(hi)
    center = enclosure.mid()
    threshold = penalty_ratio(center)[0] + (enclosure-center)*penalty_ratio(enclosure)[1]
    assert arb(".00364249896012") < threshold < arb(".00364249896016")
    assert outside_ratio < threshold < arb(1)/256
    print("Continuum excess at root minus .003", (weight(root-delta)-best).str(10),
          "local critical coefficient", critical.str(14))
    print("Global critical coefficient",threshold.str(14),"at",lo.union(hi).str(21),
          "with",accepted,"outside cells")


if __name__ == "__main__":
    verify()
