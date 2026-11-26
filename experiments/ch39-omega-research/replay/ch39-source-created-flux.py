"""Endogenous creation of degree three along the actual eight-zero clock edge network."""
from pathlib import Path
import runpy
from flint import acb,arb

fit=runpy.run_path(str(Path(__file__).with_name("ch39-hierarchy-in-zero-flux.py")))
src=fit["src"]
W,s=src["W"],src["shell"]
height=acb.zeta_zero(8).imag
p,q=[s[name]*height for name in ("a","b")]
fixed=list(map(arb,("1.180",".847","1.956","1.489")))


def degrees(k):
    loops={}
    for f in src["all_faces"]():
        value=acb(1)
        for i,j in zip(f,f[1:]+f[:1]):
            value*=-k[W.index(i),W.index(j)]
        loops[f]=[value]
    margin=min((arb.pi()-abs(v[0].arg())).lower() for v in loops.values())
    return src["degrees"](loops,0),margin


if __name__=="__main__":
    noise=arb(0,arb("1e-6"))
    model=fit["model"](p+noise,q+noise)
    _,_,original=fit["evaluate"](model,[arb(0)]*4,(0,0))
    _,_,k=fit["evaluate"](model,[x+noise for x in fixed],(0,0))
    assert degrees(original)[0]==(0,[0])
    result,margin=degrees(k)
    assert result==(3,[3]) and margin>arb(".001")
    eig=sorted([abs(z) for z in fit["symbols"](-fit["block"](k,fit["ii"],fit["oo"]))],key=lambda z:float(z.mid()))
    assert eig[2]>eig[1] and eig[1]>eig[0]
    print("created degree; angular cut margin; light/heavy, middle/heavy",
          result,margin.str(12),(eig[0]/eig[2]).str(12),(eig[1]/eig[2]).str(12))
