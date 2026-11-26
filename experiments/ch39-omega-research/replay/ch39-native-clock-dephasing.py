"""Exact dephasing by a single native clock with finite positive rational weights."""
from fractions import Fraction as Q
from itertools import combinations
from pathlib import Path
from runpy import run_path

root=Path(__file__).resolve().parent
amb=run_path(str(root/"ch39-phase-support-ambiguity.py"))
base=amb["base"]
word=base["word"]
unit=word["DEN"]


def exponents(ids):
    prefix=base["data"]("ch39-tick-separators.json")["prefix"]
    rate=base["rates"][0]
    return [-43235114615559591407865398948523*(((prefix[i-1][0]>>112)*rate+2**311)//2**312)
            for i in ids]


def cover(ids):
    powers=exponents(ids)
    gaps={abs(i-j) for i,j in combinations(powers,2)}
    assert len(gaps)==len(ids)*(len(ids)-1)//2
    phases={v:word["power"](v) for v in powers}
    width=0
    longest=0
    for i,j in combinations(powers,2):
        c,si=word["multiply"](phases[i],word["conjugate"](phases[j]))
        # Doubling a point from the open right semicircle first enters the
        # left semicircle in a finite dyadic number of turns for this roster.
        m=1
        while c[0]>0:
            m*=2
            c,si=word["multiply"]((c,si),(c,si))
            assert m<=2**16
        assert c[1]<0
        width+=2*m
        longest=max(longest,m)
    print(len(ids),"heights,",len(gaps),"distinct differences; positive rational filter has degree",
          width,"; largest doubling",longest)
    print("Bit length of sum of clock-exponent differences:",sum(gaps).bit_length())
    return width


def main():
    q,p=35937,35888
    assert q-p==49 and q%2==1
    from math import gcd
    assert gcd(2*p,q)==1
    assert all(word["power"](n)[0][0]>0 for n in range(31))
    assert word["power"](31)[0][1]<0 and word["power"](60)[1][0]>0 and word["power"](61)[1][1]<0
    roster=sorted(set(amb["upper"]["time"]["roster"]()[0]))
    assert cover(roster)==16016
    assert cover(sorted(roster+list(amb["outside"])))==33368


if __name__=="__main__":
    main()
