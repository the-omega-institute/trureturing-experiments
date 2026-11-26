"""Exact degree and number-field budget for isolated four-clock Wilson minima."""
from pathlib import Path
from math import isqrt
import runpy

m=runpy.run_path(str(Path(__file__).with_name("ch39-shell-algebraic-lock.py")))
den=m["lo"].denominator+m["hi"].denominator
# 4! Vol(unit l1 ball in dimension 4)=16. The chosen prime valuation is four.
def first_degree(field):
    n=(den+64*field-1)//(64*field)
    v=isqrt(isqrt(n))
    return v if v**4==n else v+1
for field in (1,2,64,2**20):
    degree=first_degree(field)
    assert 64*field*(degree-1)**4<den<=64*field*degree**4
    print("field degree",field,"lowest Fourier degree allowed",degree)
print("minimum coefficient field degree when total clock traversals <=5",
      (den+64*5**4-1)//(64*5**4))
