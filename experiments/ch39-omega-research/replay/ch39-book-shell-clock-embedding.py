"""Map the RH book's R=17 shell to the clock subgraph of the R=867 shell and certify independence."""
from fractions import Fraction
from itertools import product,combinations
from math import isqrt,prod
from pathlib import Path
import runpy
from flint import arb,acb

m=runpy.run_path(str(Path(__file__).with_name("ch39-minimal-triangle-readout.py")))
q,pi=m["q"],arb.pi()
v=sorted((p for p in product(range(5),repeat=3) if sum(x*x for x in p)==17),
         key=lambda p:(0 not in p,p))
digits=(1,5,13,23,29)
image=[tuple(digits[k] for k in p) for p in v]
assert len(v)==9 and image==list(m["sh"]["sites"])[:9]
assert all(sum(k*k for k in p)==867 for p in image)
def chord(p,i,j):
    return sum((x-y)**2 for x,y in zip(p[i],p[j]))
classes={}
for i,j in combinations(range(9),2):
    a,b=chord(v,i,j),chord(image,i,j)
    classes.setdefault(a,set()).add(b)
    assert (a<=14)==((i,j) in m["edges"])
    assert (a==12)==(b==724) and (a==14)==(b==804) and (a==6)==(b==244)
assert classes[2]=={32,200}
assert all(not (i<9 and j<9) for (i,j),d in m["edges"].items() if d==12)

angles=[]
classes_sqrt=[]
for R,ds in ((17,(6,12,14)),(867,(244,724,804))):
    for d in ds:
        cs=1-Fraction(d,2*R)
        assert (2*cs).denominator>1
        classes_sqrt.append(cs.numerator**2-cs.denominator**2)
        angles.append(q(cs).acos())
# The six conjugations flipping one square root each exist.
for size in range(1,7):
    for subset in combinations(classes_sqrt,size):
        n=prod(subset)
        assert n<0 or isqrt(n)**2!=n
def cross(x,y):
    return [x[(i+1)%3]*y[(i+2)%3]-x[(i+2)%3]*y[(i+1)%3] for i in range(3)]
def dot(x,y):
    return sum(a*b for a,b in zip(x,y))
n1,n2=cross(v[0],v[7]),cross(v[0],v[8])
axis=q(dot(n1,n2))/q(dot(n1,n1)*dot(n2,n2)).sqrt()
axis_sq=Fraction(dot(n1,n2)**2,dot(n1,n1)*dot(n2,n2))
assert axis>0 and axis<1
f=[angles[1]-angles[0],angles[2]-angles[0],angles[4]-angles[3],angles[5]-angles[3]]
h=acb.zeta_zero(719).imag
p,t=h*f[0],h*f[1]
s17,t17=p.sin()**2,t.sin()**2
mixed=1-(p.cos()*t.cos()+axis*p.sin()*t.sin())**2
eps=q(".001")
def probability(s):
    assert eps<s<1-eps
    return 2/pi*((s+eps).sqrt().asin()-(s-eps).sqrt().asin())
fraction=probability(s17)*probability(t17)
print("embedded squared chord classes", sorted((a,sorted(b)) for a,b in classes.items()))
print("R=17 squared cosine of the axes",axis_sq,"clock frequencies",*[x.str(12) for x in f[:2]])
print("R=17 triangular deficits at 719",*[x.str(13) for x in (s17,t17,mixed)])
print("relative density of original-shell two-window test in any positive-density new-shell clock box",
      fraction.str(13))
