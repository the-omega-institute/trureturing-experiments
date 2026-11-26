"""Uniform word-length bound from four actual triangle turns in the nine-point shell."""
from fractions import Fraction as Q
from itertools import permutations
from math import comb, prod
from pathlib import Path
from runpy import run_path
import json

DEPTH = 136
BASIS = [
    [-18321227,-7222099,18536782,-4851601,71520150],
    [-5760212,-25462894,-122016359,-34541419,-4372279],
    [-63316251,134287772,63638582,-85907930,-9185845],
    [-74391082,73467030,-74599884,126093211,46614659],
    [188469430,196739150,-70358130,-38450196,98093486],
]
TRIANGLES = [(Q(76,77),Q(-3,77)),(Q(208,217),Q(-15,217)),
             (Q(8,9),Q(-1,9)),(Q(13,21),Q(-4,21))]


def det(matrix):
    n=len(matrix)
    def term(perm):
        inversions=sum(perm[a]>perm[b] for a in range(n) for b in range(a+1,n))
        return (-1)**inversions*prod(matrix[i][j] for i,j in enumerate(perm))
    return sum(term(p) for p in permutations(range(n)))


def phase_cell(z,ops):
    # The sign of the imaginary part gives the next binary half-turn digit.
    value=0
    for _ in range(DEPTH):
        lo,hi=z[1]
        assert lo>0 or hi<0
        value=2*value+int(hi<0)
        z=ops["multiply"](z,z)
    return value


def word(exponents,ops):
    out=ops["ONE"]
    for m,g in zip(exponents,TRIANGLES):
        v=tuple(map(ops["interval"],g))
        if m<0:
            v=ops["conjugate"](v)
        k=abs(m)
        while k:
            if k%2:
                out=ops["multiply"](out,v)
            v=ops["multiply"](v,v)
            k//=2
    return out


def main():
    root=Path(__file__).resolve().parent
    ops=run_path(str(root/"ch39-shortest-triangle-clock-words.py"))
    law=run_path(str(root/"ch39-two-prime-height-law.py"))
    table=json.loads((root/"ch39-five-tick-record.json").read_text())
    rates=law["rate_numerators"]()
    codes=[phase_cell(tuple(map(ops["interval"],z)),ops) for z in TRIANGLES]
    period=1<<DEPTH
    assert all((sum(m*a for m,a in zip(row[:4],codes))-row[4])%period==0 for row in BASIS)
    determinant=det(BASIS)
    assert abs(determinant)==period
    cof=[(-1)**i*det([row[:4] for k,row in enumerate(BASIS) if k!=i]) for i in range(5)]
    assert [sum(c*r[j] for c,r in zip(cof,BASIS)) for j in range(5)]==[0,0,0,0,determinant]
    half_width=[(1+sum(abs(row[j]) for row in BASIS))//2 for j in range(5)]
    budget=sum(half_width[:4])
    assert budget==713215622 and Q(16*(1+sum(half_width)),period)<Q(121,10000*2**96)
    count=lambda length:sum(2**j*comb(4,j)*comb(length,j) for j in range(5))
    assert count(18567041)<2**96<=count(18567042)
    longest=0
    worst=Q(0)
    den=1<<table["bits"]
    unit=ops["DEN"]
    for idx,h,_,a,b in table["rows"]:
        targets=[-43235114615559591407865398948523*((v*h+(den<<23))//(den<<24)) for v in rates]
        out=[]
        for power in targets:
            code=phase_cell(ops["power"](power),ops)
            coeff=[(Q(code*c,determinant)+Q(1,2))//1 for c in cof]
            cache=[sum(m*row[j] for m,row in zip(coeff,BASIS)) for j in range(5)]
            assert abs(cache[4]-code)<=half_width[4] and all(abs(m)<=b for m,b in zip(cache[:4],half_width))
            longest=max(longest,sum(map(abs,cache[:4])))
            out.append(word(cache[:4],ops))
        u,v=out
        ca,sa=ops["standard"](a,den)
        cb,sb=ops["standard"](b,den)
        for z,t in zip([u,v,ops["multiply"](ops["conjugate"](u),v),ops["multiply"](u,v)],
                       [ca,cb,ca*cb+sa*sb,ca*cb-sa*sb]):
            worst=max(worst,max(abs(Q(x,unit)-t) for x in z[0]))
    assert worst<Q(1,2**96)
    print("Uniform budget",budget,"actual maximum for 92 words",longest,"cosine error / 2^-96",float(worst*2**96))


if __name__=="__main__":
    main()

