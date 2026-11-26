"""Explicit finite two-clock concentration at the first positive zeta height."""
from flint import arb, ctx


def main():
    ctx.prec=320
    bandwidth=arb(1)/64
    frequency_cutoff=arb(2)**68
    order=2**112
    p,q=325919355854421968365,205632218873398596256
    beat=abs(p*arb(2).log()-q*arb(3).log())
    assert beat>arb("8.90e-23") and beat<arb("8.91e-23")
    assert beat*(frequency_cutoff+15)<1
    fejer_error=4*arb.pi()/(order*beat*bandwidth)
    assert fejer_error<arb("1.75e-9")
    # Tail by averaging unitary translations up to 2/B; the transform convention
    # is the one with raw spectral measure |zeta(1/2+ix)/(1/2+ix)|² dx/(2 pi).
    B=frequency_cutoff
    s=2/B
    T=B/4
    difference2=2*(2*s*s*T+2*s*(1+T.log())+1/T)+4*s*s
    tail=80*(1+B.log())/B
    assert difference2<tail
    # On |x-gamma|<=delta/2 the zeta derivative projection exceeds .7.
    core=bandwidth**3/40000
    assert arb(".7")**2/(2*arb.pi()*226*12)>arb(1)/40000
    # Include pole cancellation (factor between 9/17 and 25/17), compact
    # smoothing and dyadic rounding of every Fejer coefficient.
    rho=fejer_error+arb(2*order)*arb(2)**(-255)
    inside=(arb(9)/17*arb(".99")*(1-rho))**2*core
    outside=(arb(25)/17)**2*(2*rho**2+(1+rho)**2*tail)
    assert outside/(inside+outside)<arb("0.00001")
    assert 2*p*order+8<2**183 and 2*q*order<2**183
    print("Finite two-clock degree",order,"gives first-zero band probability > .99999.")
    print("Upper bound on leakage", (outside/(inside+outside)).str(8))


if __name__=="__main__":
    main()
