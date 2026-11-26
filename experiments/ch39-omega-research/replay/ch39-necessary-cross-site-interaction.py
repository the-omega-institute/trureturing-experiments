"""Norm cost of leaving a three-site diagonal Yukawa texture."""
from flint import arb

q=arb(".165")
mu=q/(1-q)
up=arb(".00316")/arb(".8573")
down=arb(".000271")/arb(".01392")
budget=arb(".014")


def direction(r,epsilon):
    denom=1-2*mu-(1+mu)*(r+epsilon)
    assert denom>0
    return ((2*(1+mu)).sqrt()*mu*(r+epsilon)+(1-2*mu)*epsilon)/denom


same=max((direction(up,e)+direction(down,budget-e)).upper() for e in (arb(0),budget))
different=mu+arb(2).sqrt()*same
assert same<arb(".035") and different<arb(".98")
print("same-site and distinct-site ceilings at combined residual 0.014",
      same.str(12),different.upper().str(12))
