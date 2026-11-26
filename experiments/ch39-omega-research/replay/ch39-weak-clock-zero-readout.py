"""An ordinary neutral triangular Wilson trace labels the remaining signed angle ambiguity."""
from pathlib import Path
import runpy
from flint import arb

m=runpy.run_path(str(Path(__file__).with_name("ch39-clock-weak-cubic.py")))
q,pi=m["q"],arb.pi()
s,t=m["s"],m["t"]
p,v=m["pent"]["target"]
axis=m["pent"]["axis"]
edges={(i,j):kind for i,j,_,_,kind in m["raw"]}
assert edges[0,7]==724 and edges[0,8]==804 and edges[7,8] not in (724,804)
cp,cq=p.cos(),v.cos()
sp,sq=p.sin(),v.sin()
direct=2-(cp*cq+axis*sp*sq)**2
other=2-(cp*cq-axis*sp*sq)**2
difference=4*abs(axis)*(s*(1-s)*t*(1-t)).sqrt()
assert difference>0 and (abs(direct-other)-difference).contains(0)
u=s+t
z=s*t
trace=13+2*u
det=11*(1+u+z)
# Joint limiting density of sum of mass squares and their product.
density=1/(11*pi**2*abs(s-t)*(z*(1-u+z)).sqrt())
# Finite resolution around the unordered spectral invariants includes both
# clock labelings. A rigorous Haar measure bound follows by integrating
# the inverse-Jacobian enclosure on this rectangle entirely inside the support.
eps=q(".001")
xx=trace+arb(0,eps)
yy=det+arb(0,eps)
uu=(xx-13)/2
vv=yy/11-1-uu
disc=uu**2-4*vv
assert disc>0
l=(uu-disc.sqrt())/2
h=(uu+disc.sqrt())/2
assert l>0 and h<1
mass=4*eps**2/(11*pi**2*(disc*vv*(1-uu+vv)).sqrt())
print("normalized mixed triangular Wilson factor, opposite sign",direct.str(13),other.str(13))
print("signed-branch separation",difference.str(13))
print("sum/product density at 719",density.str(13))
print("zero density in the ±0.001 sum/product rectangle",mass.lower().str(13),mass.upper().str(13))
