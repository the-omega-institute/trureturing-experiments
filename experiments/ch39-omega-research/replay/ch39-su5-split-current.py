"""One vectorlike 10 per family suffices to remove the light-light diquark gauge vertex."""
from pathlib import Path
from itertools import combinations
import runpy
from flint import acb_mat, arb

root = Path(__file__).resolve().parent
hol = runpy.run_path(str(root/'ch39-ground-holonomy.py'))
phase = hol['phase']
charges = [-2,-2,-2,3,3]
weights = [charges[i]+charges[j] for i,j in combinations(range(5),2)]
assert sorted(set(weights))==[-4,1,6] and weights.count(1)==6
eigen = {w:phase**w for w in set(weights)}
assert min(abs(eigen[1]-eigen[w]).lower() for w in (-4,6))>arb('1.98')
def projector(value):
    return (value-eigen[-4])*(value-eigen[6])/((eigen[1]-eigen[-4])*(eigen[1]-eigen[6]))
for w in weights:
    assert abs(projector(eigen[w])-int(w==1))<arb('1e-40')
# For one chiral 10, add 10_B and its conjugate. The copy-space mass row
# is m(1-P,P); the massless column can be taken as (P,1-P).
for w in weights:
    p = int(w==1)
    assert (1-p)*p+p*(1-p)==0 and (1-p)**2+p**2==1
# The center twist of any commuting discrete charge drops out of this weight relation.
# Baryon triality charges in ninth roots of unity for Q, u^c, e^c.
assert 2*1-(-4)-6==0
assert (2*1-(-1)-0)%9==3
print('diquark-current copy overlap zero; SU(5)-orbit triality discrepancy',3,'modulo',9)
# Independently split the conjugate fundamental into the D and L branches.
d_phase,l_phase = phase**2,phase**-3
assert abs(d_phase-l_phase)>arb('1.98')
for w in (2,-3):
    p=(phase**w-l_phase)/(d_phase-l_phase)
    assert abs(p-int(w==2))<arb('1e-40')

# Restrict the standard renormalizable SU(5) vertices to the chosen light subspaces.
light = {'A':{'Q'},'B':{'u','e'},'C':{'d'},'D':{'L'}}
up_doublet = [('Q','u'),('u','Q')]
up_triplet = [('Q','Q'),('u','e'),('e','u')]
down_doublet = [('Q','d'),('e','L')]
down_triplet = [('u','d'),('Q','L')]
def surviving(first,second,channels):
    return [(a,b) for a,b in channels if a in light[first] and b in light[second]]
assert surviving('A','B',up_doublet)==[('Q','u')]
assert surviving('A','C',down_doublet)==[('Q','d')]
assert surviving('B','D',down_doublet)==[('e','L')]
assert not surviving('A','B',up_triplet)
assert not surviving('A','C',down_triplet) and not surviving('B','D',down_triplet)
for copy in light.values():
    for a,b in (('Q','u'),('Q','e'),('d','L')):
        assert not {a,b}<=copy
print('light doublet Yukawa types survive independently; light-pair vector and triplet vertices vanish')
