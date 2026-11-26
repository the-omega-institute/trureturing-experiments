"""A finite common-threshold gauge group needs a large adjoint loop parameter."""
from pathlib import Path
import runpy
from flint import arb

m = runpy.run_path(str(Path(__file__).with_name("ch39-curvature-spectrum-rigidity.py")))
q,A = m["q"],m["weyl"]
pi = arb.pi()
x = q(559)
k = 4*A/pi.sqrt()
# Every Dirichlet floor is at least its continuous argument minus one.
running = k*(x-1)/2 - q("135.5")*x.log()
inv = q("40.9906")
required_index = 3*(running-2*pi*inv)/(11*x.log())
factor = (2+required_index)/(4*pi*inv)
assert required_index > q("2081.36")
assert factor > q("4.0445")
print("extra weak adjoint index required",required_index.lower().str(14))
print("at least one simple factor's adjoint loop parameter",factor.lower().str(14))
