"""Integer exterior-character convolutions for the amplified SU(5) face and weak insertion."""
from itertools import combinations_with_replacement
from math import comb, factorial, lcm, prod
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx, fmpz_mat

full = run_path(str(Path(__file__).with_name("ch39-full-color-obstruction.py")))
q = full["q"]
DEGREE, WEAK_DEGREE = 55, 13


def alternant(shape):
    powers = [r+4-i for i,r in enumerate(shape)]
    columns = ((1,0),(1,1),(1,2),(-1,0),(-1,1))
    return int(fmpz_mat([[comb(k,j)*x**(k-j) if j<=k else 0 for x,j in columns]
                         for k in powers]).det())


def vandermonde(shape):
    return prod(shape[i]-shape[j]+j-i for i in range(5) for j in range(i+1,5))


def multiplicity(shape, vander, power, facts):
    if shape[0]>power:
        return 0
    # Dual Cauchy coefficient, then sum over determinant shifts on SU(5).
    num = vander*prod(facts[power+i] for i in range(5))
    den = prod(facts[power+i-r]*facts[r+4-i] for i,r in enumerate(shape))
    term = num//den
    assert term*den==num
    total = term
    for k in range(power-shape[0]):
        num = term*prod(power+i-r-k for i,r in enumerate(shape))
        den = prod(r+k+5-i for i,r in enumerate(shape))
        term = num//den
        assert term*den==num
        total += term
    return total


def moment_table():
    shapes = []
    dims = []
    zero = (0,)*5
    base = alternant(zero)
    for rows in combinations_with_replacement(range(WEAK_DEGREE+3),4):
        shape = (*reversed(rows),0)
        numerator = alternant(shape)
        chi = numerator//base
        assert chi*base==numerator
        if chi:
            v = vandermonde(shape)
            dims.append(v//vandermonde(zero))
            shapes.append((shape,v,chi))
    common = lcm(*dims)
    facts = [factorial(i) for i in range(DEGREE+8)]
    multiplicities, scaled = [], []
    for (shape,v,chi),dim in zip(shapes,dims):
        row = [multiplicity(shape,v,j,facts) for j in range(DEGREE+3)]
        multiplicities.append(row)
        scaled.append([a*chi*(common//dim) for a in row[:WEAK_DEGREE+3]])
    table = fmpz_mat(multiplicities).transpose()*fmpz_mat(scaled)
    moments = [[arb(int(table[i,j]))/common for j in range(table.ncols())] for i in range(table.nrows())]
    assert moments[0][0]==1 and moments[1][0]==moments[0][1]==2 and moments[2][0]==8
    assert all(table[i,j]==table[j,i] for i in range(WEAK_DEGREE+3) for j in range(i))
    return moments


def observe(moments, index, power):
    t = q(1)/15
    strength = 1/(acb.zeta_zero(index).imag**2+q(1)/4)
    alpha = strength/(1-2*(1-t)*strength)
    X, Y = 32*power*alpha, 32*power*alpha*t
    err = (X+Y).exp()*((X+Y)**(DEGREE+1)/factorial(DEGREE+1)
                       +Y**(WEAK_DEGREE+1)/factorial(WEAK_DEGREE+1))
    coeffs = [arb(1)]
    for n in range(DEGREE):
        coeffs.append(coeffs[-1]*(power-n)*alpha/(n+1))
    offsets = ((0,0),(1,0),(0,1),(2,0),(1,1),(0,2))
    sums = [arb(0) for _ in offsets]
    for n,cn in enumerate(coeffs):
        for j in range(min(n,WEAK_DEGREE)+1):
            weight = cn*comb(n,j)*(-t)**j
            for h,(a,b) in enumerate(offsets):
                sums[h] += weight*moments[n-j+a][j+b]
    boxed = [s+arb(0,32**(a+b)*err) for s,(a,b) in zip(sums,offsets)]
    A,B,AA,AB,BB = [v/boxed[0] for v in boxed[1:]]
    assert boxed[0]>1
    mean = A-t*B-2*(1-t)
    variance = AA-2*t*AB+t*t*BB-(A-t*B)**2
    print(index, "ordinary and weak-inserted traces", (A-2).str(10), (B-2).str(10),
          "face trace",mean.str(11),"variance",variance.str(10),flush=True)
    return mean,variance


def verify():
    ctx.prec = 180
    moments = moment_table()
    high = full["high"]
    power = (high["geo"]["BOOK_PLANCK"]/(1000*q(full["L0"]).exp()))**2/(8*arb.pi())
    values = [observe(moments,i,power) for i in range(11485,11491)]
    gaps = [a[0]-b[0] for a,b in zip(values,values[1:])]
    centers = (1540434462,1540136483,1539681752,1539486151,1539073109,1538878040)
    for (mean,var),v in zip(values,centers):
        assert abs(mean-q(v)/10**9) < q("0.000000002") and q("9.55")<var<q("9.57")
    assert min(gaps)>q("0.00019506")
    print("Smallest successive gap",min(gaps).str(10))


if __name__ == "__main__":
    verify()
