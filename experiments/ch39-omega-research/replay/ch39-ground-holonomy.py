"""The based U(2) holonomy closure of the eighth-zero noncollinear ground state."""
from pathlib import Path
import runpy
from flint import acb

root = Path(__file__).resolve().parent
m = runpy.run_path(str(root/'ch39-gaussian-momentum.py'))
q = m['real']
edges = {(i,j):(d,z) for i,j,d,z in m['shell']['edges']}
def walk(vertices):
    phase = acb(1)
    letters = []
    for i,j in zip(vertices,vertices[1:]):
        length,z = edges[min(i,j),max(i,j)]
        phase *= z if i<j else z.conjugate()
        if length in (724,804):
            letters.append((length,1 if i<j else -1))
    return phase,letters

assert walk((0,7,6,0))[1] == [(724,1)]
assert walk((0,8,6,0))[1] == [(804,1)]
phase,letters = walk((0,6,7,8,6,0))
assert not letters
assert m['shell']['R'] == 3*17**2
# Each squared unit phase lies in Q(sqrt(-3)); its torsion has order dividing six.
separation = abs(phase**12-1)
assert separation > 0
c = q('-.1855').union(q('-.127'))
halftrace = 1-2*(m['su']*m['sv'])**2*(1-c*c)
assert halftrace > q('-.354') and halftrace < q('-.325')
assert m['x'] > q('.545') and m['x'] < q('.546') and m['y'] > q('.15')
# Check the cosines for projective rotation order at most five.
from flint import arb
for order in range(1,6):
    for j in range(order):
        assert not m['x'].overlaps(abs((arb.pi()*j/order).cos()))
print('commutator half-trace',halftrace.lower().str(12),halftrace.upper().str(12))
print('central torsion exclusion |phase^12-1|',separation.lower().str(12))
