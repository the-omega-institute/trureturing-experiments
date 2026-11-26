"""Top and bottom phase-volume constants, and explicit small-rotation bounds for CP-coupled alignment."""
from pathlib import Path
from itertools import combinations
import runpy
from flint import arb
from fractions import Fraction as F

axes=runpy.run_path(str(Path(__file__).with_name("ch39-shell-four-family-axes.py")))
geom=runpy.run_path(str(Path(__file__).with_name("ch39-shell-family-alignment.py")))
d2,ball=axes["d2"],axes["ball"]
M=sum(d2.values())
b=[sum(c for J,c in d2.items() if i in J) for i in range(4)]
prod=F(1)
for t in b:
    prod*=t
coeff=1/(2*arb.pi()**2*ball(prod).sqrt())
minimal=ball(min(d2.values())).sqrt()
k=4*arb(int(geom["transverse"]))/geom["den"]
# If two normals rotate together, rotate their cross product as one; its norm is <=1.
# Hence each determinant of a mixed triple changes by at most dR=2 sin(theta/2).
print("maximum squared CP strength",ball(M).str(15),
      "upper-tail coefficient for squared norm",coeff.str(14))
print("smallest absolute triple determinant",minimal.str(14),
      "alignment lower penalty per squared operator-norm rotation", (k/2).str(14))
assert ball(M)>arb(".0056") and (k/2)>3
# For any lambda, the continuous envelope with maximal amplitudes beats every
# exactly coplanar contender, since the unconstrained E minimum has M>0.
