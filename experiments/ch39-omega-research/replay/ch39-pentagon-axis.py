"""First positive fundamental Wilson loops capable of selecting a noncollinear ground state."""
from pathlib import Path
from itertools import combinations, permutations
import runpy
from flint import arb, acb

root=Path(__file__).resolve().parent
m=runpy.run_path(str(root/'ch39-gaussian-momentum.py'))
q=m['real'];sh=m['shell'];pi=arb.pi()
edges={}
for i,j,ds,z in sh['edges']:
    edges[i,j]=ds,z
    edges[j,i]=ds,z.conjugate()
def cycle_info(sites):
    ds=[]
    phase=acb(1)
    signs=[]
    for i,j in zip(sites,sites[1:]+sites[:1]):
        length,z=edges[i,j]
        ds.append(length);phase*=z
        if length in (724,804): signs.append((length,1 if i<j else -1))
    return ds,phase.real,signs
# All length-four cycles with at least two instances of each clock are entirely on clock edges.
clock_vertices=sorted({i for i,j in edges if edges[i,j][0] in (724,804)})
squares=[]
for members in combinations(clock_vertices,4):
    for rest in permutations(members[1:]):
        s=(members[0],)+rest
        if s[1]>s[-1] or not all((i,j) in edges for i,j in zip(s,s[1:]+s[:1])):continue
        ds,phase,word=cycle_info(s)
        if ds.count(724)!=2 or ds.count(804)!=2:continue
        assert ds[0]==ds[2] and ds[1]==ds[3]
        assert all(word[j][1]==(-1)**j for j in range(4))
        assert phase>q('.999')
        squares.append(s)
assert len(squares)==3
ds,phase,word=cycle_info((0,7,4,6,8))
assert word==[(724,1),(804,-1),(724,1),(804,-1)]
assert phase<q('-.98')
c0=-m['x']*m['y']/(m['su']*m['sv'])
curvature=-4*phase*(m['su']*m['sv'])**2
print('pentagon center and curvature',c0.lower().str(12),c0.upper().str(12),
      curvature.lower().str(12),curvature.upper().str(12))
flow=runpy.run_path(str(root/'ch39-spin-phase-closure.py'))
t0,t1=q('26.98'),q('27.025')
T=t0.union(t1)
# Restrict the previous eighth-zero rectangle to c<=0.09 and certify the new phase face.
cos_end=q('.09')
xs=flow['spectrum'](sh['a']*m['g'],q('.00001'),True,sh['b']*m['g'],cos_end.acos())
for _ in range(2):
    time,_,_,_=flow['run'](T,xs,8192,return_target=True)
    assert T.overlaps(time)
    T=T.intersection(time)
_,margin,_,_=flow['run'](T,xs,16384)
print('restricted crossing face',margin.lower().str(12),margin.upper().str(12))
assert margin>0 and c0>cos_end
# In the continuum every crossing certified on this subrectangle can be made a unique minimum.
curv_budget=q('.15488')/(c0+q('.5'))
assert curv_budget>q('.19632')
print('combined curvature lower', (curv_budget-q('.19632')).lower().str(12))
