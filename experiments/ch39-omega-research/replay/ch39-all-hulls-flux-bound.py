"""A uniform flux obstruction for all three-dimensional convex hulls of shell vertices."""
from itertools import combinations
from pathlib import Path
import runpy
from flint import arb

model = runpy.run_path(str(Path(__file__).with_name('ch39-symmetric-shell-actions.py')))
geometry, sites = model['geometry'], model['sites']
area = model['area']
faces = geometry['boundary'](range(19))
total_area = sum(area(f) for f in faces)
assert total_area < arb('2080.244')
def face_data(vertices):
    faces = geometry['boundary'](vertices)
    areas = [area(f) for f in faces]
    angles = [geometry['lifted'](f,model['negative']) for f in faces]
    return sum(areas),sum(v*v/a for v,a in zip(angles,areas))
small = (6,9,11,13)
zero = (0,6,11,14)
a0,e0 = face_data(zero)
a1,e1 = face_data(small)
assert a0<arb('71.054') and e0<arb('.00000669')
assert a1<arb('51.979') and e1<arb('1.551')
assert geometry['degree'](geometry['boundary'](zero),model['negative'])[0]==0
# Classify areas with exactly the same integer squared cross-product signatures.
def signature(nodes):
    result=[]
    for i,j,k in combinations(nodes,3):
        u=geometry['diff'](sites[j],sites[i])
        v=geometry['diff'](sites[k],sites[i])
        n=geometry['cross'](u,v)
        result.append(geometry['dot'](n,n))
    return tuple(sorted(result))
tiny=[]
ties=[]
for nodes in combinations(range(19),4):
    p,q,r,t = [sites[i] for i in nodes]
    norm=geometry['cross'](geometry['diff'](q,p),geometry['diff'](r,p))
    if geometry['dot'](norm,geometry['diff'](t,p))==0:
        continue
    values=signature(nodes)
    if values==signature(small):
        tiny.append(nodes)
    elif values==signature(zero):
        ties.append(nodes)
    else:
        assert sum(arb(x).sqrt()/2 for x in values)>a0
assert tiny==[(6,9,11,13),(7,10,12,15),(8,16,17,18)] and len(ties)==6
for nodes in tiny:
    assert geometry['degree'](geometry['boundary'](nodes),model['negative'])[0]==-1
    assert face_data(nodes)[1]<arb('1.551')
for count,field in [(2,arb('.005')),(3,arb('.0075'))]:
    # Minkowski bounds apply without assuming the external field is constant.
    large = 2*arb.pi()*count/total_area.sqrt()-field*total_area.sqrt()
    reference = e0.sqrt()+field*a0.sqrt()
    tiny_bound = 2*arb.pi()*count/a1.sqrt()-field*a1.sqrt()
    tiny_reference = arb('1.551').sqrt()+field*a1.sqrt()
    assert large>reference and tiny_bound>tiny_reference
    print('flux, maximum field, large-area and tiny-area norm margins',
          count,field.str(5),(large-reference).lower().str(12),(tiny_bound-tiny_reference).lower().str(12))
