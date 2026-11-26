"""Construct moment-matched right scalars at a one-loop crossing with Q/M=400."""
from pathlib import Path
import runpy

from flint import arb

root=Path(__file__).resolve().parent
m=runpy.run_path(str(root/"ch39-slope-matched-one-loop.py"))
q=m["q"]
r=q(400)
tau,required=m["g"]["base"]
low=[v/r**2 for v in m["xs"]]
left=-sum(y.log() for y in low)/2
WR=left-required
yR=(-2*WR/18).exp()
assert 0<WR and 0<yR and yR<1 and -yR.log()/2<tau
lowR=18*(-yR).exp()
# The common measure has the low 18 and at least this many additional
# levels in [y19,2], by boundary-zero upper trials on the finest allowed mesh.
y19=m["edge19"]/r**2
count=sum(n for n,factor in m["factors"] if factor/(1-m["h"]**2)**2/r**2 < 2)
assert 1<y19 and y19<2 and count>18
high_count=count-18
S_lower=sum((-y).exp() for y in low)+high_count*(-q(2)).exp()-lowR
# (y-1)*exp(-y) increases on [1,2]; all higher leftovers contribute nonnegatively.
excess=-18+high_count*(y19-1)*(-y19).exp()+(1-yR)*lowR
assert S_lower>0 and excess>1/q(1).exp()
# Hypercharge gauge coupling at the projected intersection; pure fermion total exposure one.
gg=m["g"]["g"]
unified=4*arb.pi()*gg["G0"][0]-(q("41/10")*tau+q("2/15")*left+q("8/15")*WR+q("8/3"))/(2*arb.pi())
assert unified>1
for name,val in (("Q GeV",1000*tau.exp()),("left logs",left),("right logs",WR),
                 ("common low right squared mass / Q^2",yR),("right high S lower",S_lower.lower()),
                 ("right high M-S lower",excess.lower()),("inverse unified alpha",unified)):
    print(name,val.lower().str(12),val.upper().str(12))
