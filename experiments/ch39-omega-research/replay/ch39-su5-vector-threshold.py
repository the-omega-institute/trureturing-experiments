"""One-loop broken-vector threshold direction and size in the specified SU(5) flux model."""
from pathlib import Path
import runpy
from flint import arb, acb_mat, fmpz_mat

root = Path(__file__).resolve().parent
g = runpy.run_path(str(root/'ch39-su5-shell-flux.py'))
shell,q,rho = g['shell'],g['q'],g['rho']
neutral = fmpz_mat(19,19)
for i,j,ds,z in shell['edges']:
    neutral[i,i] += 1
    neutral[j,j] += 1
    neutral[i,j] = neutral[j,i] = -1
trees = fmpz_mat([[neutral[i,j] for j in range(18)] for i in range(18)]).det()
ratio = g['matrix'].det().real/(19*trees)
assert ratio>0 and 1/(2*rho)>q(74)
error = (1+108*rho).log()
log_product = ratio.log()+(-19*error).union(18*error)
r = (shell['theta'](804)/shell['theta'](12))*q('.999').union(q('1.001'))
exposure = r.log()-log_product/2
# Positive sign in alpha^-1 for b^heavy=-(7/2)T when running upwards.
delta23 = q(7)/(4*arb.pi())*exposure
print('neutral tree count',trees)
print('relative full spectral determinant at zero (dimensionless)',
      (ratio*19/g['volume']).lower().str(12),(ratio*19/g['volume']).upper().str(12))
print('sub-Q log exposure contrast',exposure.lower().str(12),exposure.upper().str(12))
print('vector threshold in alpha2^-1-alpha3^-1',delta23.lower().str(12),delta23.upper().str(12))
assert exposure>0 and delta23>0
