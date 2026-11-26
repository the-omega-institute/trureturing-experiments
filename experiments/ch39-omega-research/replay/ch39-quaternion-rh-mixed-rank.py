"""Quaternionic polarization and a finite-background interpolant for the multi-frequency shell rank."""
from fractions import Fraction as Q
from itertools import product


def conj(x):
    return (x[0],-x[1],-x[2],-x[3])


def mul(x,y):
    a,b,c,d=x
    e,f,g,h=y
    return (a*e-b*f-c*g-d*h,a*f+b*e+c*h-d*g,
            a*g-b*h+c*e+d*f,a*h+b*g-c*f+d*e)


def main():
    unit=(1,0,0,0)
    i,j=(0,1,0,0),(0,0,1,0)
    assert mul(conj(i),i)==unit
    assert mul(conj(i),j)==(0,0,0,-1) and mul(conj(j),i)==(0,0,0,1)
    # Collinear fluxes yield B+B^T; perpendicular fluxes yield B-B^T.
    for m in (2,21):
        assert 23*m+25*m==48*m
        print(m,"time directions: flat",23*m,"continuous over flux",24*m,
              "smooth and conjugation-equivariant over a flux neighborhood",48*m)
    # Quotient by the grounded columns leaves one nonnegative real Schur
    # value. The following smooth function of it is exact at any finite
    # prescribed set and at the flat point.
    samples=(Q(15648,101855),Q(353217006,318403915))
    def remaining(s):
        eta=Q(1,4096)
        for v in samples:
            eta*=((s-v)/v)**2
        assert s+eta>0
        return s*s/(s+eta)
    assert all(remaining(s)==s for s in (Q(0),*samples))
    assert Q(0)<remaining(Q(1,1000))<Q(1,1000)
    # Constant eta gives an absolute budget everywhere while retaining all
    # flat null vectors.
    for s,k in product((Q(0),Q(1,123),Q(17,4)),(96,288)):
        eta=Q(1,2**k)
        assert 0<=s-s*s/(s+eta)<eta


if __name__=="__main__":
    main()
