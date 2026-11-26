"""Cosmological coefficient needed to maintain the reference mean volume."""
from flint import arb, fmpq


def verify():
    lo, hi = 30484476, 30484477
    root5 = arb(5).sqrt()
    # Shape Laplace transform at b compared with each a in the smaller interval.
    gap = lo*(fmpq(6,250) - (fmpq(11,10)+fmpq(4,3))/(50*root5))-690
    assert gap > 50000
    assert (fmpq(5,4)-fmpq(11,10))*50000 == 7500
    assert fmpq(11,10)/(fmpq(4,3)-fmpq(5,4)) < 16
    assert fmpq(121,100)*(1-fmpq(16,2**16)) > 1
    max_per_cell = (12000/root5).sqrt()
    assert root5*(fmpq(15,lo)+fmpq(6284,1000)*max_per_cell)/2 < 600
    # With M_P^2 a^2=4 tau lambda and bare cosmological curvature K a^2=L/lambda.
    lower = 1/(50*4*fmpq(6284,1000)*hi)
    assert lower > fmpq(26,10**12)
    print("Radial exponent margin per cell above 50000:", gap)
    print("Lower bound on bare curvature in Planck units:", arb(lower))


if __name__ == "__main__":
    verify()
