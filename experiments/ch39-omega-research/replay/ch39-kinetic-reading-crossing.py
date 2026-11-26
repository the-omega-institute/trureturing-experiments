"""Rational enclosures linking a finite kinetic reading, a normalization, and a mass."""
from fractions import Fraction as Q
from math import prod
from pathlib import Path
from runpy import run_path

root=Path(__file__).resolve().parent
previous=run_path(str(root/"ch39-vote-weighted-integers.py"))
base=previous["base"]
B,rounded,multiply,absolute=[previous[k] for k in ("B","round_out","multiply","absolute")]
A=previous["A"]
ETA=Q(1,128)
TOL=Q(1,2**272)


CALIBRATIONS = (
    (3,67948801317971,19,1200,2150),
    (4,89827319625526,18,1400,2400),
)


def precision_cases():
    specs=[]
    for j,mass,bits,outer,witness in CALIBRATIONS:
        center=ETA/(j*j)
        delta=Q(outer,10**9)
        min_mass=mass-mass//10**6
        max_mass=mass+mass//10**6
        specs.extend([(center*(1+delta),min_mass,1),(center*(1-delta),max_mass,-1),
                      (center*(1-Q(witness,10**9)),max_mass+1,1)])
    return specs


def check_precision(model,rows):
    levels=base["data"]("ch39-minimum-tick-finite-inputs.json")["levels"]
    for j,mass,bits,outer,witness in CALIBRATIONS:
        center=ETA/(j*j)
        lo=Q(levels[j],2**288)
        target=(B(1)-B(lo,lo+Q(1,2**288)))/B(ETA)
        delta=Q(outer,10**9)
        tau=Q(1,2**bits)
        assert model.reading(center*(1+delta)).hi<target.lo-tau
        assert model.reading(center*(1-delta)).lo>target.hi+tau
        counter=model.reading(center*(1-Q(witness,10**9)))
        assert target.hi-2*tau<counter.lo<=counter.hi<target.lo+2*tau
        print("RH offset",j,"minimal calibration bits for one ppm:",bits)
    checked={(Q(r["c"]),r["mass"],r["sign"]) for r in rows}
    assert set(precision_cases())<=checked


class ReadingBounds:
    def __init__(self):
        data=base["data"]
        record=data("ch39-cutoff-threshold-inverses.json")[1]
        weights=base["exact_fit"](record)
        prefix=data("ch39-tick-separators.json")["prefix"]
        def ordinate(i):
            h=prefix[i-1][0]
            return rounded(B(Q(2**800,(h+1)**2),Q(2**800,h*h)))
        self.z=[ordinate(i) for i in range(1,45)]
        self.poles=[rounded(B(1)/v) for v in self.z[:5]]
        self.ref=[(ordinate(i),w) for i,w in weights.items()]
        self.p,self.modes,self.q2,self.target=base["physical_inputs"]()
        self.product=multiply(*self.z).hi
        self.tails=[]
        for a in A:
            top=multiply(*(v+B(Q(1,a)) for v in self.z))
            den=prod(abs(1-Q(k,a)) for k in A if k!=a)
            self.tails.append(rounded(top/B(den)).hi)
        for i,v in enumerate(self.z[:5]):
            other=self.z[:i]+self.z[i+1:]
            factor=multiply(*(B(1)+B(a)*v for a in A))/multiply(*(absolute(v-t) for t in other))
            sensitivity=factor.hi*sum(h/(Q(1,a)+v.lo) for a,h in zip(A,self.tails))
            assert weights[i+1].lo-2*multiply(factor,*other).hi-TOL*sensitivity>Q(".005")

    def response(self,y):
        """Return the average of 1/(1-y*z), for a point off all poles."""
        assert self.z[-1].hi*max(0,y)<1
        terms=[rounded(B(y)*v*w/(B(1)-B(y)*v)) for v,w in self.ref]
        denominator=multiply(*(absolute(B(1)-B(y)*v) for v in self.z)).lo
        mult=abs(y*prod(y+a for a in A))/denominator
        if mult:
            tail=2*self.product/min(1,1-y*self.z[-1].hi)
            tol=TOL*sum(h/abs(1+y/a) for a,h in zip(A,self.tails))
            error=mult*(tail+tol)
        else:
            # A calibrated node. This routine subtracts two calibrated tables.
            error=TOL if y else 0
        return rounded(sum(terms,B(1-error,1+error)))

    def reading(self,c):
        return rounded((B(1)-self.response(-ETA/c))/B(ETA))

    def verify_sample(self,row):
        c=Q(row["c"])
        if row.get("mass") is None:
            top=B(row["Y"]) if "Y" in row else self.poles[2]
            mass2=self.q2/(B(c)*top)
            pure=sum(v.hi<top.lo for v in self.poles) if "Y" in row else 2
        else:
            mass2=B(row["mass"]**2)
            top=self.q2/(B(c)*mass2)
            pure=sum(v.hi<top.lo for v in self.poles)
            assert pure<len(self.poles) and self.poles[pure].lo>top.hi
        assert (mass2*B(8)).lo>self.q2.hi and (mass2*B(c)*self.poles[0]).lo>1000**2
        higgs=B(125**2)/mass2
        exposure=multiply(*[(top/t)**25 for t in self.poles[:pure]])
        active=25*pure
        for lam,n,block in zip(self.modes,base["MODES"],row["roots"],strict=True):
            assert len(block)==pure+1
            def f(y):
                return B(c*y)-higgs-lam+lam*(B(1)-self.p)*self.response(y)
            for j,pair in enumerate(block):
                l,h=[Q(v,2**row.get("root_bits",75)) for v in pair]
                assert l>(self.poles[j-1].hi if j else 0) and h<self.poles[j].lo
                assert f(l).hi<0 and f(h).lo>0,(str(c),row["mass"],j,float(l),float(f(l).hi),float(f(h).lo))
                assert (B(c*l)*mass2).lo>1000**2
                if h<top.lo:
                    exposure=multiply(exposure,(top/B(l,h))**n)
                    active+=n
                else:
                    assert top.hi<l
        if row.get("mass") is None and "Y" not in row:
            assert active==119
        gap=exposure**11-self.target
        assert (gap.lo>0 if row["sign"]==1 else gap.hi<0),(str(c),row["mass"],row["sign"],float(gap.lo/self.target.lo),float(gap.hi/self.target.lo))
        assert active==row["active"]
        return active


def main():
    model=ReadingBounds()
    rows=base["data"]("ch39-kinetic-reading-roots.json")
    for row in rows:
        model.verify_sample(row)
    lo,hi=[Q(s) for s in (".00082657012849023",".00082657012849025")]
    assert [(Q(r["c"]),r["sign"]) for r in rows if r["mass"] is None]==[(lo,-1),(hi,1)]
    low=model.reading(hi).lo
    high=model.reading(lo).hi
    print("Critical kinetic reading in",float(low),float(high))
    assert Q("1.8780365378188")<low<high<Q("1.8780365378190")
    mass2=model.q2/(B(lo,hi)*model.poles[2])
    assert 69559539507206**2<mass2.lo<mass2.hi<69559539507209**2
    for j,lo_m,hi_m,active in [(3,67948801317967,67948801317975,119),(4,89827319625521,89827319625531,144)]:
        c=ETA/(j*j)
        assert [(r["mass"],r["sign"],r["active"]) for r in rows if Q(r["c"])==c]==[(lo_m,1,active),(hi_m,-1,active)]
    # On the entire interval there are exactly three open mixed roots per
    # block. Relative to the crossing, Y=Q0^2/(c M^2) lies between 624 and 635.
    c_min,c_max=Q(1,2048),Q(1,1152)
    assert model.q2.lo/(c_max*67948801317975**2)>624
    assert model.q2.hi/(c_min*89827319625521**2)<635
    assert model.poles[1].hi<624<model.poles[2].lo
    assert model.poles[2].hi<635<model.poles[3].lo
    assert c_max*635<8 and model.q2.lo*Q(1,10**5)/635>1000**2
    higgs_max=Q(125**2*635,model.q2.lo)*c_max
    for lam in model.modes:
        def secular(y,c):
            return B(c*y)-lam+lam*(B(1)-model.p)*model.response(y)
        assert secular(Q(1,10**5),c_max).hi<0 and Q(125**2*635,model.q2.lo)<Q(1,10**5)
        assert secular(Q(624),c_min).lo>higgs_max and secular(Q(635),c_max).hi<0
    check_precision(model,rows)
    print("Both on-grid mass brackets and the third-pole transition checked.")


if __name__=="__main__":
    main()
