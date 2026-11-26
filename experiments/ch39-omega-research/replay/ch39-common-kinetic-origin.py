"""Calibrate the 520-bit targets and isolate two nearby critical kinetic readings."""
import json
from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path
import sys
from flint import arb,arb_series,ctx

root=Path(__file__).resolve().parent
source=run_path(str(root/"ch39-common-kinetic-reading.py"))
base=source["base"]
ball=run_path(str(root/"ch39-fixed-poles-exposure.py"))["ball"]


def levels():
    ctx.prec,ctx.cap=1100,3
    t=arb_series([arb(1)/2,1])
    xi=t*(t-1)/2*((-t/2)*arb.pi().log()).exp()*(t/2).gamma()*t.zeta()
    mass=xi[2]/xi[0]
    level=run_path(str(root/"ch39-extra-offset-cutoffs-origin.py"))["level"]
    return [int((level(j,mass)*2**520).floor().unique_fmpz()) for j in base["OFFSETS"]]


def isolate(model,centers=False,root_bits=190,pad_bits=160):
    ctx.prec=360
    p=ball(model.p).mid()
    Y=ball(model.poles[2]).mid()
    H=125**2*Y/ball(model.q2).mid()
    target=ball(model.target).mid()
    poles=[ball(v).mid() for v in model.poles]
    support=[(ball(z).mid(),ball(w).mid()) for z,w in model.ref]
    old=base["data"]("ch39-kinetic-reading-roots.json")[0]
    def evaluate(c,save=False):
        exposure=25*sum((Y/v).log() for v in poles[:2])
        slope=arb(0)
        blocks=[]
        for lam,n,block in zip(map(ball,model.modes),base["MODES"],old["roots"]):
            values=[]
            for l,r in block:
                y=arb(l+r)/2**76
                for _ in range(5):
                    delta=lam*(1-p)*sum(w*z/(1-y*z) for z,w in support)
                    fy=c+lam*(1-p)*sum(w*z/(1-y*z)**2 for z,w in support)
                    y=(y-(c*(y-H)-lam*p+y*delta)/fy).mid()
                exposure+=n*(Y/y).log()
                slope+=n*(1-H/y)/fy
                if save:
                    step=y*arb(2)**(-pad_bits)
                    values.append([int(((y-step)*2**root_bits).floor().unique_fmpz()),int(((y+step)*2**root_bits).ceil().unique_fmpz())])
            blocks.append(values)
        return 11*exposure-target.log(),11*slope,blocks
    c=arb(".00082657012849024")
    for _ in range(4):
        v,d,_=evaluate(c)
        c=(c-v/d).mid()
    # Rounding slack keeps the two signs independent of point root guessing.
    if centers:return c,evaluate
    mid=int((c*2**165).floor().unique_fmpz())
    critical_reading=sum(w*z/(c+z/128) for z,w in support).mid()
    rows=[]
    for k,sign in ((mid-450,-1),(mid+451,1)):
        rational=Q(k,2**165)
        _,_,blocks=evaluate(arb(k)/2**165,True)
        rows.append({"c":str(rational),"mass":None,"sign":sign,"active":119,"roots":blocks,"root_bits":190})
    return rows,critical_reading


def main():
    high=root/"ch39-kinetic-reading-high-levels.json"
    new=levels()
    if "--write" in sys.argv:high.write_text(json.dumps(new)+"\n")
    else:assert json.loads(high.read_text())==new
    models=[source["FixedSpectrum"](t) for t in (Q(1,4),Q(3,4))]
    samples,critical=zip(*(isolate(model) for model in models))
    for bits in range(140):
        n=int((critical[0]*2**bits).floor().unique_fmpz())+1
        if arb(n)/2**bits<critical[1]:break
    witness={"bits":bits,"numerator":n,"samples":samples}
    path=root/"ch39-common-kinetic-witness.json"
    if "--write" in sys.argv:path.write_text(json.dumps(witness)+"\n")
    else:assert json.loads(path.read_text())==json.loads(json.dumps(witness))
    print("Critical readings",*[v.str(50) for v in critical])


if __name__=="__main__":
    main()
