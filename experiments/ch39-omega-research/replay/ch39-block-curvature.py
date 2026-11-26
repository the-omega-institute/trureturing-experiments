"""Exact combinatorics and interval bounds for additive blocks of four-simplices."""
from itertools import combinations

from flint import arb, fmpq


def verify():
    count = [6 if j-i in (1,4) else 4 for i,j in combinations(range(5),2)]
    assert count.count(4) == count.count(6) == 5
    square_deficit = sum(max(1,fmpq(d,2)-1)**2/fmpq(d) for d in count)
    assert square_deficit == fmpq(55,12)
    assert 24*square_deficit*60 == 6600
    root5 = arb(5).sqrt()
    fraction = root5/(4*400**2*6600)
    assert fraction > fmpq(1,2**31)
    lo,hi = 30484476,30484477
    correction = 15*root5/(2*lo)
    trace = (fmpq(1,50)-correction)/(4*fmpq(6284,1000)*hi)
    assert trace > fmpq(26,10**12)
    print("Minimum product of anomalous block fractions:",fraction)
    print("Volume-weighted signed trace curvature lower bound in Planck units:",trace)


if __name__ == "__main__":
    verify()
