"""A lower bound on the parallel-transport gauge-boson mass form of the explicit shell background."""
from collections import deque
from pathlib import Path
import runpy

from flint import acb_mat

root = Path(__file__).resolve().parent
color = runpy.run_path(str(root/'ch39-eighth-color-holonomy.py'))
m = runpy.run_path(str(root/'ch39-gaussian-momentum.py'))
q = m['real']
basis = color['traceless_basis']()
unit = color['identity'](8)
gram = acb_mat(8,8)
for g in color['generators']:
    images = [g*b*color['adjoint'](g) for b in basis]
    diff = acb_mat([[color['hilbert_inner'](a,b) for b in images] for a in basis])-unit
    gram += color['adjoint'](diff)*diff

def positive_definite(matrix):
    # Sylvester's criterion for the exact Hermitian matrix enclosed by these entries.
    for k in range(1,matrix.nrows()+1):
        block = acb_mat([[matrix[i,j] for j in range(k)] for i in range(k)])
        if not block.det().real > 0:
            return False
    return True

lower = q(0)
for j in range(1,401):
    trial = q(j)/100
    if not positive_definite(gram-unit*trial):
        break
    lower = trial
assert lower>0
weak = 4*min((m['su']**2).lower(),(m['sv']**2).lower())*(1-q('.1855'))
links = [[] for _ in range(19)]
for i,j,_,_ in m['shell']['edges']:
    links[i].append(j)
    links[j].append(i)
distances = {0:0}
queue = deque([0])
while queue:
    i = queue.popleft()
    for j in links[i]:
        if j not in distances:
            distances[j] = distances[i]+1
            queue.append(j)
assert len(distances)==19
radius_sum = sum(distances.values())
rho = m['shell']['rho']*m['shell']['theta'](804)
def gap(value):
    endpoint = q(37)*6/value + 2*radius_sum
    with_edges = (1+2*max(map(len,links))*rho)*endpoint+2*rho
    return 1/with_edges
print('color and weak loop bounds',lower.str(12),weak.lower().str(12))
print('graph distance sum',radius_sum)
print('dimensionless squared-mass gap, color and weak',gap(lower).lower().str(12),gap(weak).lower().str(12))
assert gap(lower)>q('.001') and gap(weak)>q('.005')
