"""A minimal-rank special-unitary flux background with the global standard-model centralizer."""
from pathlib import Path
import runpy
from flint import acb_mat

root = Path(__file__).resolve().parent
m = runpy.run_path(str(root/'ch39-gaussian-momentum.py'))
q = m['real']
shell = m['shell']
matrix = acb_mat(19,19)
for i,j,ds,z in shell['edges']:
    matrix[i,i] += 1
    matrix[j,j] += 1
    matrix[i,j] = -z**5
    matrix[j,i] = matrix[i,j].conjugate()
def positive_shift(s):
    # Each Hermitian leading principal minor must be positive.
    for n in range(1,20):
        block = acb_mat([[matrix[i,j] - (s if i==j else 0) for j in range(n)] for i in range(n)])
        if not block.det().real>0:
            return False
    return True
lo,hi = 0,18000
while hi-lo>1:
    mid = (lo+hi)//2
    if positive_shift(q(mid)/1000):
        lo = mid
    else:
        hi = mid
# A negative leading principal minor at the upper endpoint certifies that the shift is not PSD.
upper = q(hi)/1000
det = (matrix-acb_mat([[upper*int(i==j) for j in range(19)] for i in range(19)])).det().real
assert det<0
lower = q(lo)/1000
# For a general field, the endpoint energy is bounded above by its full edge energy.
rho = shell['rho']*shell['theta'](804)
mass_floor = 1/((1+36*rho)/lower+2*rho)
green = matrix.inv()
diagonal = [green[i,i].real for i in range(19)]
print('static charged ground squared mass',lower.str(12),upper.str(12))
print('all edge models lower charged squared mass',mass_floor.lower().str(12))
print('endpoint current Green diagonals',min(v.lower() for v in diagonal).str(12),max(v.upper() for v in diagonal).str(12))
# The common normalized four-dimensional zero mode is constant on the connected graph.
# The interior mass quadrature on each unit edge lies between zero and its full edge mass.
aggregate = sum(shell['rho']*shell['theta'](ds) for i,j,ds,z in shell['edges'])
volume = q(19).union(q(19)+aggregate)
ratio = shell['theta'](804)/shell['theta'](12)*q('.999').union(q('1.001'))
mass = 1000*q('26.994').union(q('27.008')).exp()/ratio
site_scale = mass/(volume*min(v.lower() for v in diagonal).union(max(v.upper() for v in diagonal))).sqrt()
print('conditional same-site current exchange scale GeV',site_scale.lower().str(12),site_scale.upper().str(12))
assert mass_floor>q('1.3724') and all(q('.08113')<v<q('.13703') for v in diagonal)
