"""Prime support, baseband beat and direct RH-zero windows for the shared-reference loss floor."""
from pathlib import Path
from runpy import run_path
from flint import arb,acb

home=Path(__file__).resolve().parent
geom=run_path(str(home/"ch39-shell-third-shell.py"))
weil=run_path(str(home/"ch39-shell-local-prime.py"))
a,B=weil["a"],weil["B"]
pairs=[((17,0),(29,1)),((17,1),(26,1)),((26,0),(29,0))]
freq=[(geom["shells"][r][2][i],geom["shells"][s][2][j]) for (r,i),(s,j) in pairs]
bands=[x for p,q in freq for x in (4*p,4*q,2*abs(p-q),2*(p+q))]
small=[v for v in bands if v<B]
assert len(small)==1
gap=small[0]
assert gap>2*a and all(v+B<arb(5).log() for v in bands)
assert all(abs(v-arb(k).log())>B for v in bands for k in (2,3,4))
amp=weil["density"](gap,False)/(8*weil["f0"])
print("detuning",gap.str(16),"asymptotic amplitude",amp.str(16),"period",(2*arb.pi()/gap).str(15))

A=arb(3200)
count=A.zeta_nzeros().unique_fmpz()
assert count is not None
targ=(25,255,719)
ts={n:acb.zeta_zero(n).imag for n in targ}
totals={n:[arb(0),arb(0)] for n in targ}
s=A+arb(1)/2
xi=1/s+1/(s-1)-arb.pi().log()/2+(s/2).digamma()/2
for start in range(1,int(count)+1,96):
    for z in acb.zeta_zeros(start,min(96,int(count)-start+1)):
        h=z.imag
        floor=sum(((p*h).sin()**2-(q*h).sin()**2)**2 for p,q in freq)/2
        for n,T in ts.items():
            w=(a*(h-T)).sinc()**4+(a*(h+T)).sinc()**4
            totals[n][0]+=w
            totals[n][1]+=w*floor
for n,T in ts.items():
    D,N=totals[n]
    tail=(2*A*xi/(a**4*(A-T)**4)).upper()
    added=arb(0).union(tail)
    den=weil["denominator"](T)
    assert (D+added).contains(den)
    local=(N+3*added/2)/den
    approx=arb(3)/8-amp*(gap*T).cos()
    print("window center",n,"certified floor mean",local.lower().str(12),local.upper().str(12),
          "leading beat",approx.str(13))
