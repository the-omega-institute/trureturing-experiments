"""Effect of a dynamical R=17 shell on the already selected R=867 zero."""
from pathlib import Path
import runpy
from flint import arb

m=runpy.run_path(str(Path(__file__).with_name("ch39-book-shell-clock-embedding.py")))
q,pi=m["q"],arb.pi()
c=m["axis"]
L=m["s17"]+m["t17"]+m["mixed"]
p,t=m["p"],m["t"]
sp,cp=p.sin_cos()
st,ct=t.sin_cos()
z=cp*ct+c*sp*st
gp=2*sp*cp-2*z*(-sp*ct+c*cp*st)
gt=2*st*ct-2*z*(-cp*st+c*sp*ct)
assert L>0 and L<1 and gp*gp+gt*gt>0

def excess(p):
    s=p.sin()**2
    A=(1-s+(1-c*c)**2*s*s/4).sqrt()
    B=1+(3-c*c)*s/2-L
    if B>A:
        return q(0)
    assert B<A and B>0
    return ((A*A-B*B).sqrt()-B*(B/A).acos())/pi
def check():
    # At fixed p the q-average is decreasing on p in [0,pi/2].
    n=8192
    terms=[excess(pi*j/(2*n)) for j in range(n+1)]
    mean=(sum(terms[1:]).lower()/n).union(sum(terms[:-1]).upper()/n)
    coefficient=(q(57429)*mean.lower()).lower().union((q(58472)*mean.upper()).upper())
    assert mean>0
    print("R=17 triangle energy and gradient at 719",L.str(13),gp.str(13),gt.str(13))
    print("mean positive part of the triangle energy deficit",mean.lower().str(12),mean.upper().str(12))
    print("density slope of lower-energy zeros per positive coupling",coefficient.lower().str(12),coefficient.upper().str(12))
if __name__=="__main__":
    check()
