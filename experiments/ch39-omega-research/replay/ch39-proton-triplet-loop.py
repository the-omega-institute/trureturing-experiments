"""Finite Euclidean tests of radiative mixing after orthogonal tree matching."""
from fractions import Fraction as Q
from math import comb
from pathlib import Path
from runpy import run_path

root = Path(__file__).resolve().parent
old = run_path(str(root/'ch39-eight-time-pole-control.py'))
B, node, rounded = (old[k] for k in ('B','node','rounded'))
mass_ratio = node(1)/node(2)


def integrated_log(q, order=58):
    # q in [0,4]. Integrate log(1+q*t*(1-t))-log(3/2) exactly up to the
    # order shown, with |(1-2*q*t*(1-t))/3| <= 1/3 on the unit interval.
    assert 0<=q<=4
    q_grid = 2**96
    mid = Q(q*q_grid//1,q_grid)
    # Integral of u**j is (j!)**2/(2j+1)!, u=t*(1-t).
    moments = [(-2*mid)**j/((2*j+1)*comb(2*j,j)) for j in range(order+1)]
    total = Q(0)
    for n in range(1,order+1):
        power = sum(comb(n,j)*moments[j] for j in range(n+1))
        total -= power/(n*3**n)
    remainder = Q(3,2*(order+1)*3**(order+1))
    # Integrated derivative is between zero and 1/6.
    return rounded(B(total-remainder,total+remainder+(q-mid)/6))


def subtracted(x):
    if x==0: return B(0)
    low = integrated_log(x/mass_ratio.hi)
    high = integrated_log(x/mass_ratio.lo)
    return rounded(B(low.lo,high.hi)-integrated_log(Q(x)))


def main():
    k1,k4 = subtracted(1),subtracted(4)
    curvature = k4-B(4)*k1
    positive = curvature/B(3)
    minimax = curvature/B(44)
    second = (B(1)-B(1)/mass_ratio**2)/B(60)
    for name,val in [('R',mass_ratio),('curvature',curvature),
                     ('local after matching at 1 and 4',positive),
                     ('minimax over three samples',minimax),
                     ('quadratic coefficient after matching at 0',second)]:
        print(name,float(val.lo),float(val.hi))
    for value,lower in ((curvature,'.08445665630617'),(positive,'.02815221876872'),
                        (minimax,'.00191946946150'),(second,'.01326024177060')):
        assert Q(lower)<value.lo<=value.hi<Q(lower)+Q(1,10**14)
    assert 3-4+1==0 and 3*0-4*1+4==0
    assert 3*1**2+4*2**2+5**2==44
    # Saturate the lower bound on the maximum of the three full amplitudes.
    a = minimax
    b = B(-4)*minimax-k1-a
    residual = (k4+a+B(4)*b)/B(25)-minimax
    assert residual.lo<=0<=residual.hi


if __name__=='__main__':
    main()
