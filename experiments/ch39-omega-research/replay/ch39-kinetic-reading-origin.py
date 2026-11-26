"""Generate finite secular root candidates; the independent replay checks their signs."""
import json
from pathlib import Path
from runpy import run_path
from fractions import Fraction as Q
from flint import arb,ctx

root=Path(__file__).resolve().parent
s=run_path(str(root/"ch39-kinetic-reading-crossing.py"))
base=s["base"]
ball=run_path(str(root/"ch39-fixed-poles-exposure.py"))["ball"]


def candidate_rows(model,specs):
    ctx.prec=320
    support=[(ball(z).mid(),ball(w).mid()) for z,w in model.ref]
    poles=[ball(v).mid() for v in model.poles]
    q2,p=ball(model.q2).mid(),ball(model.p).mid()
    rows=[]
    for c,mass,sign,*fixed_y in specs:
        c=Q(c)
        scale=arb(c.numerator)/c.denominator
        top=arb(fixed_y[0]) if fixed_y else poles[2] if mass is None else q2/(scale*mass**2)
        mass2=q2/(scale*top)
        pure=2 if mass is None and not fixed_y else sum(v<top for v in poles)
        active=25*pure
        blocks=[]
        for lam,n in zip(model.modes,base["MODES"]):
            l=ball(lam).mid()
            block=[]
            for a,b in zip([arb(0)]+poles,poles[:pure+1]):
                for _ in range(104):
                    y=(a+b)/2
                    f=scale*y-125**2/mass2-l+l*(1-p)*(1+y*sum(w*z/(1-y*z) for z,w in support))
                    if f<0:a=y
                    else:b=y
                pad=(a+b)/2*arb("2e-17")
                block.append([int(((a-pad)*2**75).floor().unique_fmpz()),int(((b+pad)*2**75).ceil().unique_fmpz())])
                if b<top:active+=n
            blocks.append(block)
        row={"c":str(c),"mass":mass,"sign":sign,"active":active,"roots":blocks}
        if fixed_y:row["Y"]=fixed_y[0]
        rows.append(row)
    return rows


def specs():
    return [(".00082657012849023",None,-1),(".00082657012849025",None,1),
            ("1/1152",67948801317967,1),("1/1152",67948801317975,-1),
            ("1/2048",89827319625521,1),("1/2048",89827319625531,-1)]+s["precision_cases"]()


if __name__=="__main__":
    import sys
    records=candidate_rows(s["ReadingBounds"](),specs())
    path=root/"ch39-kinetic-reading-roots.json"
    if "--write" in sys.argv:path.write_text(json.dumps(records)+"\n")
    else:assert records==json.loads(path.read_text())
