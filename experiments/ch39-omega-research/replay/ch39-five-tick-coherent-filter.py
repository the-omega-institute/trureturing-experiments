"""One real two-clock wavefront for the 42 RH levels and the five-tick nine-point shell."""
from fractions import Fraction
from pathlib import Path
from runpy import run_path

from flint import acb, acb_series, arb, ctx

data=run_path(str(Path(__file__).with_name("ch39-five-tick-record-origin.py")))


def main():
    _,heights,weights=data["model"]()
    two=arb(2)
    delta,guard=two**-300,two**-285
    N,B=2**3200,two**1500
    with ctx.workprec(4096):
        n,e=(arb(3).log()/two.log()).mid().man_exp()
        ratio=(Fraction(int(n))*Fraction(2)**int(e)).limit_denominator(2**1536)
        p,q=ratio.numerator,ratio.denominator
        beat=abs(p*two.log()-q*arb(3).log())
        assert two**-1545<beat and beat<two**-1505
        assert beat*(B+1800)<1 and max(2*p*N+8,2*q*N)<2**4740
    assert all(g-h>1 for g,h in zip(heights[1:],heights))
    spill=256*arb.pi()/(N*beat*delta*guard)+2*N*two**-4500
    assert spill<two**-1060
    pole=lambda g:1-arb(8)/17*(arb(16).log()*g).cos()
    with ctx.workprec(384):
        for g,w in zip(heights,weights):
            band=g+arb(0,2*delta)
            slope=abs(acb_series([acb(arb(1)/2,g),1]).zeta()[1])
            second=2*abs(acb_series([acb(arb(1)/2,band),1]).zeta()[2])
            assert slope>arb(".75") and second<20 and band<1800
            gain=(w*(g*g+arb(1)/4)).sqrt()/(pole(g)*slope)
            assert gain>arb(".5") and gain<255
            modulation=1+arb(0,20*delta/slope)
            density=(modulation*pole(band)/pole(g))**2*(arb(1)/4+g*g)/(arb(1)/4+band**2)
            assert abs(density-1)<two**-290
    # Each signed pair's ideal mass is w*S. Near each outer edge the
    # positive Fejer convolution is at most that band's gain plus spill.
    S=2*delta**3/(3*arb.pi())
    inner=(1-guard)**3*(1-two**-290)*(1-2*spill)**2*(1-1800*two**-325)**2
    outer=(1+guard)**3*(1+two**-290)*(1+2*spill)**2
    assert inner>1-two**-282 and outer<1+two**-282
    tail=80*(1+B.log())/B
    s,T=2/B,B/4
    assert 2*(2*s*s*T+2*s*(1+T.log())+1/T)+4*s*s<tail
    outside=(arb(25)/17)**2*(2*spill**2+257**2*tail)/S
    assert outside<two**-560
    mass_error=two**-282+outside
    eps=(2*mass_error+(1+guard)*delta*outer)/(1-mass_error)
    assert eps<two**-280
    phase=5*arb(6).log()*eps
    assert 8*phase/(1-4*phase)<two**-273
    print("46 paired bands, one coherent real wavefront, order 2^3200 and positive prime exponents below 2^4740.")
    print("All 42 squared resolvent readings within 2^-280, all normalized shell inner products within 2^-273.")


if __name__=="__main__":
    main()
