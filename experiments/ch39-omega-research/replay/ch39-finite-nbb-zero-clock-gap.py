"""A dyadic dilation readout and the cost of resolving the first zero in the raw NBB spectrum."""
from flint import acb, acb_series, arb, arb_series, ctx


def main():
    ctx.prec, ctx.cap = 192,3
    s = arb_series([arb(1)/2,1])
    xi = s*(s-1)/2*((-s/2)*arb.pi().log()).exp()*(s/2).gamma()*s.zeta()
    total = xi[2]/xi[0]
    ordinates = [z.imag for z in acb.zeta_zeros(1,80)]
    weights = [1/(total*g**2) for g in ordinates]
    lag = 41*arb(2).log()
    partial = sum(p*(lag*g).cos() for p,g in zip(weights,ordinates))
    tail = 1-sum(weights)
    nbb = (2+lag)/arb(2)**(arb(41)/2)
    assert tail>0 and partial-tail-nbb>arb(".265")
    assert nbb<arb(".000021") and weights[0]/2>arb(".108")
    g=ordinates[0]
    width=arb(1)/64
    for j in range(32):
        left=g-width+2*width*j/32
        right=g-width+2*width*(j+1)/32
        stripe=acb(arb(1)/2,left.union(right))
        assert abs(acb_series([stripe,1]).zeta()[1])<1
    assert g-width>14 and 2*arb.pi()>6 and 1764*arb(".1")>13**2
    print("Dyadic lag 41: zero-clock excess over raw NBB correlation above .265.")
    print("For half-width <=1/64, putting probability .1 at the first positive zero requires normalized multiplier gain >13/width^(3/2).")


if __name__=="__main__":
    main()
