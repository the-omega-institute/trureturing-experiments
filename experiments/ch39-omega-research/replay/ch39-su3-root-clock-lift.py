"""Lift the clocked shell parallel transports to the Y=2 root doublet of a single SU(3)."""
from collections import deque
from pathlib import Path
import runpy
from flint import acb,acb_mat,arb

short = runpy.run_path(str(Path(__file__).with_name("ch39-short-arc-vacuum.py")))
m,sh,q = short["m"],short["shell"],short["q"]
source = m["w"]["source"]
axis = source["gauge"]["spectrum"]
p,t = source["clocks"]
sp,cp = p.sin_cos()
st,ct = t.sin_cos()
c,v = axis["cos_axis"],axis["sin_axis"]
U = acb_mat([[acb(cp,sp),0],[0,acb(cp,-sp)]])
V = acb_mat([[acb(ct,c*st),v*st],[-v*st,acb(ct,-c*st)]])
I = acb_mat([[1,0],[0,1]])
# A spanning tree without either clock makes both chord holonomies visible up to U(1).
adj = [[] for _ in sh["sites"]]
for i,j,d,z in sh["edges"]:
    if d not in (724,804):
        adj[i].append(j)
        adj[j].append(i)
seen = {0}
todo = deque(seen)
while todo:
    for j in adj[todo.popleft()]:
        if j not in seen:
            seen.add(j)
            todo.append(j)
assert len(seen)==19
gap = 4*(v*sp*st)**2
assert gap>0
# On an edge the SU(2) clock has determinant one. The scalar z is lifted by
# diag(w I_2, w^-2), with w^3=z. For negative real z choose angle +pi.
for i,j,d,z in sh["edges"]:
    if z.imag.is_zero() and z.real<0:
        theta = arb.pi()
    else:
        theta = z.arg()
    w = acb(0,theta/3).exp()
    clock = U if d==724 else V if d==804 else I
    lift = acb_mat(3,3)
    for a in range(2):
        for b in range(2):
            lift[a,b]=w*clock[a,b]
    lift[2,2]=1/(w*w)
    assert (lift.det()-1).contains(0)
    assert all((lift[a,b]/lift[2,2]-z*clock[a,b]).contains(0) for a in range(2) for b in range(2))
print("lifted edges",len(sh["edges"]),"ordinary tree vertices",len(seen))
print("squared norm of the clock commutator U V - V U",gap.lower().str(12),gap.upper().str(12))

# Canonical tree matching at the one-loop intersection if all 18 complex vectors
# are distinct four-dimensional gauge generators and the threshold shifts are neglected.
# alpha_U^{-1}<41 at both intersections, g_Y^2 = 12*pi/(5 alpha_U^{-1}).
gmin2=12*arb.pi()/(5*41)
loop=144*gmin2/(16*arb.pi()**2)
assert loop>q(".167")
assert 18>arb.pi()/gmin2
print("lower bound on max adjoint loop factor",loop.lower().str(12))
