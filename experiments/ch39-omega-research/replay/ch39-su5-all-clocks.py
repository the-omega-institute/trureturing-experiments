"""Arbitrary-clock obstruction in the minimally completed SU(5) Wilson model."""
from pathlib import Path
import runpy
from flint import acb, fmpz_mat

root = Path(__file__).resolve().parent
local = runpy.run_path(str(root/'ch39-su5-complete-matter.py'))
f,q = local['flow'],local['q']
edges = {(i,j):(d,z) for i,j,d,z in f['shell']['edges'] if j<9}
lap = fmpz_mat(9,9)
for i,j in edges:
    lap[i,i] += 1
    lap[j,j] += 1
    lap[i,j] = lap[j,i] = -1
trees = fmpz_mat([[lap[i,j] for j in range(1,9)] for i in range(1,9)]).det()
assert trees == 91125
walk = (0,2,8,3,1,6,4,5,7,0)
assert sorted(walk[:-1])==list(range(9))
phase = acb(1)
letters = []
for i,j in zip(walk,walk[1:]):
    ds,z = edges[min(i,j),max(i,j)]
    phase *= z if i<j else z.conjugate()
    if ds in (724,804):
        letter = (ds,1 if i<j else -1)
        if letters and letters[-1]==(ds,-letter[1]):
            letters.pop()
        else:
            letters.append(letter)
assert not letters
upper = sum(q(int(lap[j,j])).log() for j in range(9))
floors = {}
shape_lower = {}
for power in (2,3,7):
    loop = abs(1-phase**power)**2/9
    gap = loop/9
    assert gap>f['high_error']
    factor = (1-f['high_error']/gap)/((1+f['eta'])*(1+f['data']['beta']*36))
    floors[power] = gap*factor
    # The Schur complement at site zero dominates loop*I_2. The rooted
    # Dirichlet determinant dominates the square of the scalar spanning-tree count.
    shape_lower[power] = (q(int(trees))*loop).log()+9*factor.log()
    print('scalar exponent, all-mesh floor and log shape lower',power,
          floors[power].lower().str(12),shape_lower[power].lower().str(12))
base = local['dot'](local['flux_projection'],[local['D']*x/2 for x in local['g']['G0']])
flux = base
for charge,row in local['indices'].items():
    weight = q('7/3')*local['dot'](local['flux_projection'],list(map(q,row)))
    flux += weight*(shape_lower[abs(charge+1)] if weight>0 else upper)
print('Hadamard log shape upper and required vector contrast lower',upper.upper().str(12),flux.lower().str(12))
required_time = max((local['vectors']['r']/floor.sqrt()).log().upper() for floor in floors.values())
print('time sufficient for all low matter masses to lie above input',required_time.str(12))
assert required_time < q(5) and flux>q(40) and local['vectors']['exposure']<q('2.12')
