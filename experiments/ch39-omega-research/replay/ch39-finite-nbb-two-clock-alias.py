"""A uniform dyadic alias obstruction for physical zero heights."""
from flint import arb, acb, acb_series, ctx


def main():
    ctx.prec,ctx.cap=160,2
    gamma=acb.zeta_zero(1).imag
    width=arb(1)/64
    period=2*arb.pi()/arb(2).log()
    for j in range(32):
        left=gamma-width+2*width*j/32
        right=gamma-width+2*width*(j+1)/32
        x=left.union(right)
        z=acb(arb(1)/2,x)
        derivative=acb_series([z,1]).zeta()[1]
        assert abs(derivative)<1 and derivative.real>arb(".7")
        assert abs(acb(arb(1)/2,x-period).zeta())>arb(".73")
        assert x>14 and x<15 and abs(x-period)<arb("5.1")
    assert (arb(1)/4+arb("5.1")**2)/(196*arb(".73")**2)<arb(".26")
    assert arb(".26")/64**2<arb(".000064")
    # Two-clock filters may impose both real pole zeros with a factor
    # that is bounded away from zero along the unitary spectral axis.
    m=arb(4)
    a=m/(1+m*m)
    assert 1-2*a>arb(".52")
    print("Dyadic alias bound for fundamental first-zero window: probability < .26 delta^2.")
    print("The first zero is simple with Re zeta' > .7 throughout the certified interval.")


if __name__=="__main__":
    main()
