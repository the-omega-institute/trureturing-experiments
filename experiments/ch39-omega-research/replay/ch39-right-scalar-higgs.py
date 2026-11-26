"""Electroweak matching including the complementary right scalar sector at its two-loop crossing."""
from pathlib import Path
import runpy

from flint import arb

root=Path(__file__).resolve().parent
m=runpy.run_path(str(root/"ch39-right-scalar-rge.py"))
q=m["q"]
T=q("26.92").union(q("26.97"))
s=q("2.55").union(q("2.605"))
Q=1000*T.exp()
D=m["D"]
g2=4*arb.pi()/q("33.1367").union(q("33.2625"))
src=(root/"ch39-ew.py").read_text().split("interval=(q(1000)")[0]
light={}
exec(compile(src,"light_SM","exec"),light)
ev=light["evolve"]
dt=(q(1000)/q("91.1876")).log()
MZ=ev(light["low_input"],-dt)
lamb=q(".5611")/4
te=ev(MZ+[lamb,(lamb*q("248.404")**2).log()],dt)
low=ev(light["low_input"]+te[12:],q("26.92"))
target=ev(low,q(".05"),True)
ll,msq=target[12],target[13].exp()
lh=q(".02335")
f=(2*arb.pi()**2*(lh-ll)/(3*m["WL"])).sqrt()
# The complete Dirac family at r=400 is at Q/e; all further fermions are heavier.
fermion=q("36/5")*g2**2*13*(-q(2)).exp()
scalar=arb(0)
for x in m["data"]["xs"]:
    y=x/400**2
    l=y.log()
    scalar+=y*(q("78/25")*(3*l**2/2+2*l+q("7/2"))*g2**2
               +q("144/5")*(3-l)*g2*f+24*((l-2)**2+1)*f**2)
yr=(-2*s).exp()
lr=-2*s
# 2 C_Y(H) T_Y(R) = 2*(3/20)*(8/5) = 12/25.
scalar+=18*yr*q("12/25")*(3*lr**2/2+2*lr+q("7/2"))*g2**2
response=8*18*yr*(1-lr)
b=(fermion+scalar+msq*(D/Q)**2)/(D*response)
assert b>0
margin=4*lh*q(".06")-f**2
cap=max(6*lh,2*109*q(".06"),2*37*q(".005"),2*145*q(".005"))
contact=cap+(216*f**2+72*(3*f)**2+288*b**2).sqrt()
sens=1+(fermion+scalar)*(Q/D)**2/msq
quartic_ratio=(b/f)**2*18*s/(3*m["WL"])
assert margin>0 and contact<8*arb.pi()
for name,z in (("low quartic",ll),("low mass GeV^2",msq),("left f",f),("right b",b),
               ("stability",margin),("contact",contact),("relative mass sensitivity",sens),
               ("single omitted quartic term relative",quartic_ratio),("Q GeV",Q)):
    print(name,z.lower().str(12),z.upper().str(12))
