"""Exact finite-ring implementation of a bounded quaternion comparison window."""
from fractions import Fraction as Q
from itertools import permutations, product
from math import gcd


def times(x,y):
    a,b,c,d = x
    u,v,w,s = y
    return (a*u-b*v-c*w-d*s,a*v+b*u+c*s-d*w,
            a*w-b*s+c*u+d*v,a*s+b*w-c*v+d*u)


def bar(x):
    a,b,c,d = x
    return a,-b,-c,-d


def main():
    points = list(permutations((0,1,3,4)))
    a,b = points[:2]
    sigma = tuple(Q(k,26) for k in times(a,bar(b)))
    unit = (1,0,0,0)
    assert times(sigma,bar(sigma)) == unit and sigma[0]==Q(25,26)
    power = unit
    prev,current = 2,25
    for n in range(1,33):
        power = times(power,sigma)
        assert 2*power[0] == Q(current,13**n)
        assert current%13 == (-1)**n%13 and (2*power[0]).denominator == 13**n
        prev,current = current,25*current-169*prev
    # A residue comparison window for all 24 generators and inverses, including
    # a composite modulus. No primality condition enters the decoding proof.
    generators = [tuple(Q(k,26) for k in times(p,bar(a))) for p in points]
    assert all(times(x,bar(x)) == unit for x in generators)
    alphabet = sorted(set(generators+[bar(x) for x in generators]))
    for horizon in (1,2,5):
        D = 26**horizon
        modulus = 2*D+1
        assert gcd(26,modulus)==1
        reduce = lambda x:tuple(int(v)%modulus for v in x)
        def encode(x):
            return tuple(v.numerator*pow(v.denominator,-1,modulus)%modulus for v in x)
        def decode(x):
            return tuple(Q((D*v+D)%modulus-D,D) for v in x)
        if horizon<=2:
            words = product(alphabet,repeat=horizon)
        else:
            words = ((x,bar(y),x,y,bar(x)) for x in alphabet for y in alphabet)
        for word in words:
            exact = residue = unit
            for x in word:
                exact = times(exact,x)
                residue = reduce(times(residue,encode(x)))
            assert decode(residue)==exact and reduce(times(residue,bar(residue)))==unit
    print("Shortest rotor:",sigma,"; 2 Re(sigma^n) has denominator 13^n.")
    print("All shell words of lengths one and two and length-five comparison samples decode exactly.")


if __name__=="__main__":
    main()
