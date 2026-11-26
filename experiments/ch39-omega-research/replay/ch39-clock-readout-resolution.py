"""Finite-resolution number density with a labeled line centroid and a mixed triangular trace."""
from pathlib import Path
import runpy
from flint import arb

m=runpy.run_path(str(Path(__file__).with_name("ch39-weak-clock-zero-readout.py")))
q,pi=m["q"],m["pi"]
s,t=m["s"],m["t"]
eps=q("1e-7")
xx=m["trace"]+arb(0,eps)
yy=m["det"]+arb(0,eps)
u=(xx-13)/2
v=yy/11-1-u
disc=u**2-4*v
hi=(u+disc.sqrt())/2
lo=(u-disc.sqrt())/2
assert 0<lo<hi<1
move=(q(".21555")+q(".19955"))*max(abs(hi-s).upper(),abs(lo-t).upper())
centroid=q("1e-5")
assert move<centroid<q("2.3630e-5")-move
axis=m["axis"]
mixed_mean=2-(1-hi)*(1-lo)-axis**2*hi*lo
branch=2*abs(axis)*(hi*(1-hi)*lo*(1-lo)).sqrt()
near,far=mixed_mean+branch,mixed_mean-branch
assert abs(near-m["direct"])<q(".01") and abs(far-m["direct"])>q(".01")
# Exactly two of the eight equal-Jacobian sheets survive. At most a
# measure-zero algebraic set loses a cubic root's oscillator strength.
mass=eps**2/(11*pi**2*(disc*v*(1-u+v)).sqrt())
print("maximum centroid drift on each labeling",move.upper().str(13))
print("zero density in the labeled resolution box",mass.lower().str(13),mass.upper().str(13))
