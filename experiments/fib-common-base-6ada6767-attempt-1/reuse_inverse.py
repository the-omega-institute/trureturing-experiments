# Publication repair (2026-10-23): resolve original suppliers within this checkout.
"""Reproduce the existing IX/X normalized carrier and pseudo-remainder.
Reused verbatim by executing the exact symbolic portion of the private certificate.
"""
from pathlib import Path
import json
SUPPLIER=(Path(__file__).resolve().parents[1]/'fib-self-calibration/fixed-b-and-negative-child/negative-child-exhaustion.py')
namespace={}
source=SUPPLIER.read_text()
exec(compile(source[:source.index('pc=Fr(-11)')],str(SUPPLIER),'exec'),namespace)
out={}
for name in ['R','alpha','zeta']:
    P=namespace[name]
    out[name]=[[list(e),int(c)] for e,c in P.terms()]
(Path(__file__).resolve().parent/'inverse_coefficients.json').write_text(json.dumps(out,separators=(',',':'))+'\n')
print('EXACT_SYMBOLIC_REUSE_PASS')
