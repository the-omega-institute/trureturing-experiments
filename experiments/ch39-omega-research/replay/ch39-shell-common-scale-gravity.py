"""Upper bound on gravitational strength before fitting when alignment also supplies family vector masses."""
from pathlib import Path
import runpy

n=runpy.run_path(str(Path(__file__).with_name("ch39-shell-alignment-newton.py")))
q,pi,mp,Q,A,B,L,baseline=[n[t] for t in ("q","pi","mp","Q","A","B","L","baseline")]
for label,min_ratio in (("equal gF,gU",q(1)),("unitary adjoint budget",q(9)/495)):
    vmax=q("40.9906")/(9*min_ratio)
    moment=A*L*vmax*vmax.sqrt()+19*B*vmax
    fraction=(Q.upper()/mp)**2*(baseline.upper()+moment/(60*pi))
    print(label,"nu upper",vmax.upper().str(11),
          "upper fraction of measured inverse Newton constant",fraction.upper().str(13))
    assert fraction<q(".0004")
