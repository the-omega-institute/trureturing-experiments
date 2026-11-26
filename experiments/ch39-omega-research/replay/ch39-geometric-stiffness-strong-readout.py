"""以零点时钟与壳探测动量固定三角刚度，留出强耦合作读数。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

vacuum = run_path(str(Path(__file__).with_name("ch39-ten-common-triangle-vacuum.py")))
soft = vacuum["soft"]
crossing = vacuum["phase"]["crossing"]
smooth = crossing["smooth"]


def verify():
    ctx.prec = 256
    h = acb.zeta_zero(11431).imag
    geo = soft["geo"]
    r = geo["theta"](804)/geo["theta"](12)
    K = (h*h+soft["point"]["study"]["offset"])/r**2
    assert K > 1688200 and K < 1688201
    op,_,vnorm,squared,w = vacuum["setup"](h)
    x,displacement = vacuum["light"](K,op,vnorm,squared,w)
    # The same global-minimum proof applies on K >= 1688200.
    assert (arb(3).sqrt()*arb.pi()*(arb("28.5")/(2*K)).sqrt()).cos() > arb(".999")
    masses = crossing["shell_spectra"]((0,5))
    motions = arb(0,(arb(2).sqrt()*displacement).upper())
    mz,em,weak = arb("91.1876"),arb("127.930"),arb(".23129")
    pi, beta = arb.pi(), smooth["beta"]
    # Both signed effective charge-five towers share the unperturbed spectrum.
    neutral,charged = (sum(smooth["sheet"](v,r*r,True) for v in masses[q]) for q in (0,5))
    d0 = sum(smooth["sheet"](v+motions,r*r,False) for v in masses[0])
    d0 += smooth["sheet"](x.lower(),r*r,False).union(smooth["sheet"](x.upper(),r*r,False))
    assert all(v+motions > 1 for row in masses.values() for v in row)
    d5 = sum(smooth["sheet"](v+motions,r*r,False) for v in masses[5])
    # Each signed exposure lies in the same interval d5; summing nonnegative
    # index coefficients preserves this enclosure without equating the spectra.
    intercept = [a+7*(n*neutral+k*charged)/(4*pi)-2*(zw*d0+fw*d5)/(3*pi)
                 for a,n,k,zw,fw in zip((3*(1-weak)*em/5,weak*em,0),(0,2,3),(5,3,2),
                                       (arb(".1"),arb("1.5"),1),(arb("1.4"),0,arb(".5")))]
    def ref_error(mu):
        ref2 = (r*mz/mu)**2
        return arb(0,(ref2/8*(19*(35/(4*pi)+1/pi)+1/(pi*x.lower()))).upper())
    wide_err = ref_error(arb("1e13"))
    dy = beta[0]-beta[1]
    L = (intercept[0]-intercept[1]+2*wide_err)*2*pi/dy
    assert mz*L.exp() > arb("1e15")
    err = ref_error(arb("1e15"))
    L = (intercept[0]-intercept[1]+2*err)*2*pi/dy
    assert mz*L.exp() < arb("1e16")
    # d_L of a reference threshold is bounded by twice its small-momentum
    # log bound. On both endpoint scales the electroweak difference changes sign.
    assert 4*wide_err.upper() < dy/(2*pi)
    for mu,sgn in ((arb("1e15"),1),(arb("1e16"),-1)):
        assert sgn*(intercept[0]-intercept[1]-dy*(mu/mz).log()/(2*pi)+2*err)>0
    U = intercept[0]-beta[0]*L/(2*pi)+err
    strong = 1/(U-intercept[2]+beta[2]*L/(2*pi)+err)
    delta = strong-arb(".1180")
    assert arb(".117967") < strong < arb(".117973")
    assert arb("6.343e15") < mz*L.exp() < arb("6.345e15")
    assert arb("53.9049") < U < arb("53.9051")
    print("Geometric stiffness, weak-hypercharge crossing, common inverse, required alpha_s and offset",
          K.str(12),(mz*L.exp()).str(12),U.str(10),strong.str(10),delta.str(9))


if __name__ == "__main__":
    verify()
