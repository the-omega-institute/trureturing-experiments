"""Noncommuting ground state with unit pentagon cost and an actual eighth-zero crossing."""
from pathlib import Path
import runpy

root = Path(__file__).resolve().parent
flow = runpy.run_path(str(root/'ch39-spin-phase-closure.py'))
unit = runpy.run_path(str(root/'ch39-unit-pentagon.py'))
m,q = unit['model'], unit['q']
u,v = (m['shell'][s]*m['g'] for s in ('a','b'))
r = flow['r']
lo,hi = r*q('.999'), r*q('1.001')
t0,t1 = q('26.994'),q('27.008')
T = t0.union(t1)
worst = [q(1)]*4
for j in range(8):
    left = (q('-.1567')+q(j)/10000).acos()
    right = (q('-.1567')+q(j+1)/10000).acos()
    xs = flow['spectrum'](u,(left-right)/2+q('.00001'),True,v,(left+right)/2)
    faces = ((t0,lo.union(hi),0,1),(t1,lo.union(hi),0,-1),(T,lo,1,1),(T,hi,1,-1))
    for i,(t,scale,k,sign) in enumerate(faces):
        trial = t
        if k==1:
            for _ in range(2):
                target = flow['run'](trial,xs,8192,scale,return_target=True)[0]
                assert trial.overlaps(target)
                trial = trial.intersection(target)
        margin = sign * flow['run'](trial,xs,16384 if k else 4096,scale)[k]
        worst[i] = min(worst[i],margin.lower())
        assert margin>0, (j,i,margin)
print('eighth zero unit-coupling scale crossing margins',*[v.str(12) for v in worst],flush=True)
