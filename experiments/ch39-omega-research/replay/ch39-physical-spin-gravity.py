"""Ordinary boson/fermion signs, massive Proca fields and local curvature coefficients."""
from fractions import Fraction as F
from pathlib import Path
import runpy

m = runpy.run_path(str(Path(__file__).with_name("ch39-short-arc-vacuum.py")))
q = m["q"]
from flint import arb

# Scalar-flat a4 tuples are (Riemann^2, Ricci^2), modulo total derivatives.
scalar = (F(1,180),-F(1,180))
spin = (-2*(scalar[0]-F(1,96)), -2*scalar[1])
proca = (3*scalar[0]-F(1,12), 3*scalar[1]+F(1,2))
assert spin == (F(7,720),F(1,90))
assert proca == (-F(1,15),F(29,60))
# 16 signed fermionic components = eight massive two-component spinors.
solutions = [(s,v) for s in range(17) for v in range(6)
             if s+3*v==16 and s*scalar[0]+v*proca[0]+8*spin[0]==0]
assert solutions==[(10,2)]
s,v = solutions[0]
xi = (s*F(1,6)-v*F(1,2)+16*F(1,12))/s
assert xi == F(1,5)
assert s*(F(1,6)-xi)-v*F(1,2) == -F(4,3)
assert s*scalar[1]+v*proca[1]+8*spin[1] == 1
# The sign of G^{-1} is positive for bosons on the right and fermions on the left.
Qmin = 1000*m["w"]["b"]["T0"].exp()
Qmax = 1000*m["w"]["b"]["T1"].exp()
Mpl = q("1.220890e19")
lower = 2*q(".05")/(3*arb.pi())*(Qmin/Mpl)**2
upper = 2*q("1.99779")/(3*arb.pi())*(Qmax/Mpl)**2
weyl = m["r"]*sum(m["total"])/arb.pi().sqrt()
assert q("1e-11")<lower<upper<q("1e-9")
print("scalar/vector/two-component spinor counts",s,v,8,"required scalar xi",xi)
print("fraction of observed inverse G",lower.lower().str(12),upper.upper().str(12))
print("dimensionless common heat leading coefficient",weyl.str(12))
print("R=0 Ricci-squared coefficient, slope per Lambda/Q",(weyl/(16*arb.pi()**2)).str(12))
