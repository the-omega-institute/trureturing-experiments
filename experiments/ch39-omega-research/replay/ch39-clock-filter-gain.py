"""Necessary gain and direction bounds for contact plus one filtered winning clock edge."""
from pathlib import Path
import runpy
from flint import arb

work=Path(__file__).resolve().parent
strength=runpy.run_path(str(work/"ch39-clock-edge-strength-window.py"))
src=strength["src"]
k,v,g=(strength[key] for key in ("k","v","g"))
e,a,b=src["fourier_columns"](k,v)
t=strength["isolate"](strength["polys"][6,14,7],".001715",".001717")
r=max(abs(x).upper() for x in e)
amplitude=lambda x:sum(abs(z)**2 for z in x).sqrt()
# Select the heaviest Fourier eigenvalue by strict separation, then enclose the ratio.
assert all(abs(e[0])>abs(z) for z in e[1:])
delta=t*amplitude(a)*amplitude(b)/abs(e[0])
overlap=abs(sum(x.conjugate()*y for x,y in zip(a,b)))/(amplitude(a)*amplitude(b))
assert delta<arb(".0004")
inputs=[(arb(".00316")/arb(".8573"),arb(".00000629")/arb(".8573"),arb(".00000629")/arb(".00306")),
        (arb(".000271")/arb(".01392"),arb(".0000137")/arb(".01392"),arb(".0000137")/arb(".000265"))]
print("maximal relative bare update and bare direction overlap",delta.str(12),overlap.str(12))
for ratio,light,over_middle in inputs:
    for cond in (arb(1),arb("1.809")):
        gain=(1-cond*light)/(cond*ratio*delta)
        align=cond*(ratio+over_middle)/(1-cond*light)
        print("endpoint condition, required gain product, allowed direction overlap",
              cond.str(5),gain.lower().str(12),align.upper().str(12))
        assert gain>arb("50000") and align<overlap
