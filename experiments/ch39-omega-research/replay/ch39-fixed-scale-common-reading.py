"""Equal normalizer, one identical finite readout and 50 identical positive sites."""
from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path

root=Path(__file__).resolve().parent
same=run_path(str(root/"ch39-common-kinetic-reading.py"))
cross,base,B,rounded=same["cross"],same["base"],same["B"],same["rounded"]
multiply=cross["multiply"]
eta=cross["ETA"]


class ThreeFits:
    def __init__(self):
        cells=base["data"]("ch39-kinetic-reading-high-levels.json")
        self.rhs=[(v,v+int(i>0)) for i,v in enumerate(cells)]
        self.fits=[]
        for rec in base["data"]("ch39-cutoff-threshold-inverses.json"):
            mat,_=base["observations"](rec["support"])
            self.fits.append(base["exact_fit"](rec,system=(mat,self.rhs+[(0,0)]*4),target_bits=520,enclose=rounded))
        self.ids=sorted(set.union(*(set(w) for w in self.fits)))
        assert len(self.ids)==50
        prefix=base["data"]("ch39-tick-separators.json")["prefix"]
        gauge=run_path(str(root/"ch39-kinetic-reading-gauge.py"))
        pi2=gauge["pi_cell"]()**2
        d=[B(8)*pi2/B(Q(g)**2) for g in (".467453",".63811","1.0583")]
        def exp(x,bits):
            out=gauge["exp_cell"](x/B(2**bits))
            for _ in range(bits):out=rounded(out**2)
            return out
        self.q2=rounded(B(10**6)*exp((B(5)*d[0]-B(3)*d[1]-B(2)*d[2])/B(22),6))
        self.target=exp(B(109)*(d[1]-d[2])-B(Q(115,2))*(d[0]-d[1]),12)
        _,_,q2,target=base["physical_inputs"]()
        assert q2.lo<self.q2.lo<=self.q2.hi<q2.hi and target.lo<self.target.lo<=self.target.hi<target.hi
        self.poles={}
        for i in self.ids:
            h=prefix[i-1][0]
            self.poles[i]=rounded(B(Q(h*h,2**800),Q((h+1)**2,2**800)))

    def difference(self,k,y):
        # (G_k-G_2)(y)/y. Both endpoint tables solve every calibrated row.
        a,b=self.fits[k],self.fits[2]
        ids=sorted(set(a)|set(b))
        assert len(ids)==49
        nodes=(45,81,127,210,473)
        p=self.poles
        parts=[]
        for i in nodes:
            v=p[i]
            numer=multiply(a[i]-b[i],*(B(1)-v/p[j] for j in ids if j!=i))
            den=multiply(*(v+B(j*j) for j in base["OFFSETS"]),*(v-p[j] for j in nodes if j!=i))
            parts.append(multiply(numer/den,*(B(y)-p[j] for j in nodes if j!=i)))
        poly=sum(parts,B(0))
        return rounded(multiply(poly,*(B(y+j*j) for j in base["OFFSETS"][1:]))/
                       multiply(*(B(1)-B(y)/p[i] for i in ids)))

    def models(self,c):
        k=rounded(self.difference(0,-eta/c)/self.difference(1,-eta/c))
        assert max(abs(k.lo),abs(k.hi))<20
        models=[]
        for sign in (-1,1):
            shift=B(Q(sign,96))
            a=B(Q(1,3))+shift
            b=B(Q(1,3))-shift*k
            d=B(1)-a-b
            assert min(a.lo,b.lo,d.lo)>0
            w={i:rounded(sum((coef*fit.get(i,B(0)) for coef,fit in zip((a,b,d),self.fits)),B(0)))
               for i in self.ids}
            m=same["FixedSpectrum"].__new__(same["FixedSpectrum"])
            m.initialize(w,self.rhs)
            m.q2,m.target=self.q2,self.target
            models.append(m)
        return k,models



def generate(tri):
    from flint import arb
    from json import dumps
    origin=run_path(str(root/"ch39-common-kinetic-origin.py"))
    trial=Q(".000826570128490238")
    _,models=tri.models(trial)
    centers=[origin["isolate"](m,True)[0] for m in models]
    for bits in range(180):
        num=int((centers[1]*2**bits).floor().unique_fmpz())+1
        if arb(num)/2**bits<centers[0]:break
    c=Q(num,2**bits)
    k,models=tri.models(c)
    rows=[]
    for model,sign in zip(models,(-1,1)):
        _,evaluate=origin["isolate"](model,True,210,170)
        residual,_,roots=evaluate(arb(c.numerator)/c.denominator,True)
        assert sign*residual>0
        rows.append({"c":str(c),"mass":None,"sign":sign,"active":119,"root_bits":210,"roots":roots})
    (root/"ch39-fixed-scale-key.json").write_text(dumps(rows)+"\n")


def main():
    tri=ThreeFits()
    samples=base["data"]("ch39-fixed-scale-key.json")
    c=Q(samples[0]["c"])
    k,models=tri.models(c)
    assert Q(".00082657012849023")<c<Q(".00082657012849025")
    for model,row,sign in zip(models,samples,(-1,1),strict=True):
        assert Q(row["c"])==c and row["mass"] is None and row["sign"]==sign
        model.verify_sample(row)
    readings=[m.reading(c) for m in models]
    assert max(v.lo for v in readings)<=min(v.hi for v in readings)
    assert all(Q("1.87803653781886")<v.lo<=v.hi<Q("1.87803653781888") for v in readings)
    print("Shared normalizer",c,"shared readout",float(readings[0].lo))
    print("50 identical positive heights, 44 identical time rows, 4 identical phases; open counts 144 and 119.")


if __name__=="__main__":
    import sys
    if "--write" in sys.argv:generate(ThreeFits())
    main()
