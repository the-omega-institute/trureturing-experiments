"""Identify the first scalar cutoff giving all 3-by-3 Yukawa matrices and bound inverse cost."""
from pathlib import Path
from runpy import run_path
import numpy as np
from flint import arb,acb,acb_mat


m=run_path(str(Path(__file__).with_name("ch39-shell-gw-higgs-low-band.py")))
ctx=m["ctx"];ctx.prec=384
sing,family,phases,g=m["transported_spectra"]()
count=m["model"]["model"]["count"]
hs=m["hs"];adj=m["adj"]
yuk=m["yuk"]


def rows(a):
    return [[a[i,j] for j in range(a.ncols())] for i in range(a.nrows())]


def frame(k,cut):
    raw=rows(k)
    n=count(raw,cut)
    vals,V=np.linalg.eigh(np.array([[complex(v) for v in row] for row in raw]))
    pred=arb(str(vals[n]))
    lower=(cut+pred)/2
    assert pred>lower>cut and count(raw,lower)==n
    cols=[]
    for i in range(n):
        y=acb_mat([[acb(float(t.real),float(t.imag))] for t in V[:,i]])
        e=arb(float(vals[i]))
        residual=k*y-e*y
        err=hs(residual,residual).real.sqrt()/(lower-e)
        assert err<arb("1e-10")
        r=arb(0,err.upper());cols.append(acb_mat([[y[j,0]+acb(r,r)] for j in range(len(raw))]))
    return cols


if __name__=="__main__":
    lo,hi=arb("1.1046581044"),arb("1.1046581046")
    assert count(rows(sing),lo)==count(rows(sing),hi)==1
    assert count(rows(family),lo)==7 and count(rows(family),hi)==8
    sources=[(v,[yuk["identity"](3)/arb(3).sqrt()]) for v in frame(sing,hi)]
    sources.extend((v,g) for v in frame(family,hi))
    left,right=yuk["zero_basis"](1),yuk["zero_basis"](-1)
    cols=[]
    for v,G in sources:
        y=acb_mat(3,3)
        for i in range(6):
            h=acb_mat(3,3)
            for j,gen in enumerate(G):h+=v[len(G)*i+j,0]*gen
            a=acb_mat([[left[3*i+t,j] for j in range(3)] for t in range(3)])
            b=acb_mat([[right[3*i+t,j] for j in range(3)] for t in range(3)])
            y+=phases[i]*a.transpose()*h.transpose()*b
        cols.append(y)
    gram=acb_mat([[hs(a,b) for b in cols] for a in cols])
    assert count(rows(gram),arb(".0319")**2)==0
    # The numerical frame was projected into the rigorous spectral band. Its synthesis
    # map has squared norm no larger than the upper bound on the scalar-space Gram.
    norm_gram=acb_mat([[hs(a,b) if len(G)==len(J) else 0 for b,J in sources] for a,G in sources])
    assert count(rows(norm_gram),arb("1.000001"))==9
    print("first universal cutoff in",lo.str(13),hi.str(13),
          "one singlet plus eight adjoints; smallest gain >", (arb(".0319")/arb("1.000001").sqrt()).lower().str(10))
