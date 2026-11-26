"""Two exactly calibrated 49-height models at one identical dyadic kinetic reading."""
from fractions import Fraction as Q
from math import isqrt
from pathlib import Path
from runpy import run_path

root=Path(__file__).resolve().parent
cross=run_path(str(root/"ch39-kinetic-reading-crossing.py"))
base,B,rounded=cross["base"],cross["B"],cross["rounded"]
DEPTH=520


class FixedSpectrum(cross["ReadingBounds"]):
    def __init__(self,t):
        cells=base["data"]("ch39-kinetic-reading-high-levels.json")
        levels=[(v,v+int(i>0)) for i,v in enumerate(cells)]
        assert len(levels)==44 and levels[0]==(2**DEPTH,2**DEPTH)
        records=base["data"]("ch39-cutoff-threshold-inverses.json")
        fits=[]
        for rec in (records[1],records[2]):
            matrix,previous=base["observations"](rec["support"])
            for (lo,hi),(a,b) in zip(levels,previous):
                assert a*2**(DEPTH-288)<=lo<=hi<=b*2**(DEPTH-288)
            fits.append(base["exact_fit"](rec,system=(matrix,levels+[(0,0)]*4),target_bits=DEPTH,enclose=rounded))
        weights={i:B(t)*fits[1].get(i,B(0))+B(1-t)*fits[0].get(i,B(0))
                 for i in set(fits[0])|set(fits[1])}
        assert len(weights)==49
        self.initialize(weights,levels)

    def initialize(self,weights,levels):
        ids=sorted(weights)
        prefix=base["data"]("ch39-tick-separators.json")["prefix"]
        self.ref=[]
        for i in ids:
            h=prefix[i-1][0]
            z=rounded(B(Q(2**800,(h+1)**2),Q(2**800,h*h)))
            w=weights[i]
            assert w.lo>0
            self.ref.append((z,w))
        self.poles=[rounded(B(1)/z) for z,w in self.ref[:5]]
        coarse,modes,self.q2,self.target=base["physical_inputs"]()
        a,b=[B(Q(l,2**DEPTH),Q(h,2**DEPTH)) for l,h in levels[1:3]]
        ab=a+B(2)*b
        self.p=rounded((B(1)+B(Q(2,3))*ab-a-b-ab**2/(B(9)*a*b))/(B(1)-b))
        assert coarse.lo<self.p.lo<=self.p.hi<coarse.hi
        self.modes=[B(6),B(4)]
        for k in (5,3):
            lower=isqrt(k<<440)
            v=B(Q(lower,2**220),Q(lower+1,2**220))
            self.modes.extend([B(3)-v,B(3)+v])
        assert all(v.lo<=x.lo<=x.hi<=v.hi for v,x in zip(modes,self.modes))

    def response(self,y):
        return rounded(B(1)+B(y)*sum((rounded(w*z/(B(1)-B(y)*z)) for z,w in self.ref),B(0)))


def main():
    data=base["data"]("ch39-common-kinetic-witness.json")
    enclosures=[]
    for t,samples in zip((Q(1,4),Q(3,4)),data["samples"],strict=True):
        model=FixedSpectrum(t)
        for row,sign in zip(samples,(-1,1),strict=True):
            assert row["mass"] is None and row["sign"]==sign and row["active"]==119
            model.verify_sample(row)
        lo,hi=[Q(row["c"]) for row in samples]
        assert Q(".00082657012849023")<lo<hi<Q(".00082657012849025")
        enclosures.append(B(model.reading(hi).lo,model.reading(lo).hi))
    lower,upper=enclosures
    delta=upper-lower
    assert Q("1.30e-41")<delta.lo<=delta.hi<Q("1.33e-41")
    den=2**data["bits"]
    common=Q(data["numerator"],den)
    assert lower.hi<common<upper.lo
    assert (lower.lo*(den//2))//1==(upper.hi*(den//2))//1
    print("Least common dyadic denominator exponent for this pair:",data["bits"])
    print("Dyadic numerator:",data["numerator"],"open doublets at the same reading: 144 and 119.")


if __name__=="__main__":
    main()
