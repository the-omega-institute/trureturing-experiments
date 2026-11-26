"""Two eigenvalue lifts on the seven shell vertices with a fixed traceless face source."""
from itertools import product
from pathlib import Path
import runpy
from flint import arb

hull = runpy.run_path(str(Path(__file__).with_name('ch39-symmetric-shell-actions.py')))
geometry = hull['geometry']
negative = hull['negative']
area = hull['area']
R = geometry['s']['R']
b = 1/(2*arb(R))
tau = 2*arb.pi()
def branch(face):
    angle = geometry['lifted'](face,negative)
    a = area(face)
    shifts = []
    margins = []
    for target in (b*a,-b*a):
        candidates = [(abs(angle+n*tau-target),n) for n in (-1,0,1)]
        # The interval for (target-angle)/2pi guarantees these three shifts suffice.
        assert abs(target-angle)<tau
        distance,n = min(candidates,key=lambda v:float(v[0].mid()))
        assert distance<arb.pi()
        shifts.append(n)
        margins.append((arb.pi()-distance).lower())
    return shifts, min(margins)

def check(faces):
    flux,_ = geometry['degree'](faces,negative)
    choices = [branch(f) for f in faces]
    assert flux==3
    assert sorted(n for (n,k),margin in choices)==[0]*len(choices)
    assert sorted(k for (n,k),margin in choices)==[-1]*3+[0]*(len(choices)-3)
    assert 2*flux+sum(sum(n) for n,_ in choices)==3
    return min(s for n,s in choices)
big = geometry['big']
all_margins = [check(big)]
quadrilaterals = [i for i,f in enumerate(big) if len(f)==4]
for choices in product((0,1),repeat=3):
    faces = []
    for i,f in enumerate(big):
        if len(f)==4:
            if choices[quadrilaterals.index(i)]:
                f=f[1:]+f[:1]
            faces += [(f[0],f[1],f[2]),(f[0],f[2],f[3])]
        else:
            faces.append(f)
    all_margins.append(check(faces))
limits = []
for f in big:
    angle=geometry['lifted'](f,negative)
    undecorated=geometry['lifted'](f,set())
    is_odd=sum(tuple(sorted((i,j))) in negative for i,j in zip(f,f[1:]+f[:1]))%2
    if is_odd:
        assert undecorated<0
        assert (angle-undecorated-arb.pi()).contains(0)
        limits.append(((arb.pi()-angle)/area(f),len(f)))
assert len(limits)==6
lower=max(v.upper() for v,n in limits if n==4)
upper=min(v.lower() for v,n in limits if n==3)
assert lower<b<upper
print('interval of source strengths for polygon faces',lower.str(15),upper.str(15))
print('Berry source, minimum distance to branch wall',b.str(15),min(all_margins).str(15))
