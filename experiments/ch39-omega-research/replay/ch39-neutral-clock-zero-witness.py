"""A positive zeta zero with lower neutral clock energy than the fitted height."""
from pathlib import Path
import runpy
from flint import acb

m = runpy.run_path(str(Path(__file__).with_name("ch39-neutral-clock-density.py")))
q = m["q"]
n = 54542
height = acb.zeta_zero(n).imag
theta = m["sh"]["theta"]
sigma = sum((height*(theta(d)-theta(244))).sin()**2 for d in (724,804))
assert sigma<q(".00000548311")
# Global bound and the independent 719th-zero triangle lower bound.
upper = q("9.9583e-8")*sigma
fitted_lower = q("3.6473e-12")
assert 6*upper<fitted_lower
print(n,height.str(16),"normalized energy bound",upper.upper().str(12))
