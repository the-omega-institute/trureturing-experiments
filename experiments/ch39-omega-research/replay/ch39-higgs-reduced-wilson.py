"""Parallel fundamental Higgs, reduced SU(4) character spectrum, and Wilson singlet."""
from functools import cache
from itertools import product
from math import comb, factorial, prod
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

full = run_path(str(Path(__file__).with_name("ch39-full-color-obstruction.py")))
q = full["q"]
origin = (0,0,0,0)
fund = (1,0,0,0)
anti = (1,1,1,0)
six = (1,1,0,0)
degree = 74
facts = [factorial(k) for k in range(degree+6)]


def vander(shape):
    return prod(shape[i]-shape[j]+j-i for i in range(4) for j in range(i+1,4))


@cache
def multiplicities(row):
    v = vander(row)
    coeffs = []
    for n in range(degree+1):
        if n<row[0]:
            coeffs.append(0)
            continue
        num = v*prod(facts[n+i] for i in range(4))
        den = prod(facts[n+i-r]*facts[r+3-i] for i,r in enumerate(row))
        term = num//den
        total = term
        for k in range(n-row[0]):
            term = term*prod(n+i-r-k for i,r in enumerate(row))//prod(r+k+4-i for i,r in enumerate(row))
            total += term
        coeffs.append(total)
    return coeffs


def expansion(power,alpha):
    coefficients = [arb(1)]
    for j in range(degree):
        coefficients.append(coefficients[-1]*(power-j)*alpha/(j+1))
    x = 16*power*alpha
    return coefficients,x.exp()*x**(degree+1)/factorial(degree+1)


def coefficient(row,series):
    coefficients,err = series
    return sum(a*j for a,j in zip(coefficients,multiplicities(row)))+arb(0,err)


def verify():
    ctx.prec = 180
    high = full["high"]
    power = (high["geo"]["BOOK_PLANCK"]/(1000*q(full["L0"]).exp()))**2/(8*arb.pi())
    for i in range(11485,11491):
        s = 1/(acb.zeta_zero(i).imag**2+q(1)/4)
        t = q(1)/15
        alpha = 2*s/(1-2*(1-t)*s)
        series = expansion(power,alpha)
        den = coefficient(origin,series)
        l2 = coefficient(origin,expansion(2*power,alpha))/den**2
        def value(row):
            return coefficient(row,series)*vander(origin)/(vander(row)*den)
        c,b = value(fund),value(six)
        other = l2.sqrt()/100
        for dynkin in product(range(7),repeat=3):
            row = tuple(sum(dynkin[k:]) for k in range(4))
            if row not in (origin,fund,anti,six) and vander(row)<100*vander(origin):
                other=max(other,abs(value(row)))
        print(i,"fundamental",c.str(10),"six",b.str(10),"other",other.str(10),"D2",l2.str(9),flush=True)
        assert 0<other<b<c and value(anti).overlaps(c)
        remainder = l2-1-32*c**2-36*b**2
        assert q(".20017")<c<q(".20037") and q(".14350")<b<q(".14366") and other<q(".04652")
        for n,bound in ((4,".572"),(5,".033")):
            width=32*c**n+36*b**n+other**(n-2)*remainder
            variance=32*c**(2*n)+36*b**(2*n)+other**(2*n-2)*remainder
            assert width<1
            slab=14*((1+width)/(1-width))**q("6.5")*variance.sqrt()/(1-width)
            assert slab<1 and slab/(1-slab)<q(bound)
            print("faces",n,"SU(4) center area factor",(slab/(1-slab)).str(6))


if __name__ == "__main__":
    verify()
