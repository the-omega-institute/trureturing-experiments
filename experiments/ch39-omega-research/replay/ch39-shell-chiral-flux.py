"""Principal flux on the convex hull of the coordinate-reflection walls of the 19-point shell."""
from itertools import combinations, product
from pathlib import Path
import runpy
from flint import acb, arb

s=runpy.run_path(str(Path(__file__).with_name('ch39-arc-relative-determinant.py')))
pts=s['sites']
edges={(i,j):(d,z) for i,j,d,z in s['edges']}
wall=tuple(i for i,v in enumerate(pts) if len(set(v))<3)
pole=pts.index((17,17,17))
inner=tuple(i for i in wall if sum(pts[i])>=49)
def diff(a,b):
    return tuple(x-y for x,y in zip(a,b))
def cross(a,b):
    return tuple(a[(j+1)%3]*b[(j+2)%3]-b[(j+1)%3]*a[(j+2)%3] for j in range(3))
def dot(a,b):
    return sum(x*y for x,y in zip(a,b))

def boundary(nodes):
    facets={}
    for i,j,k in combinations(nodes,3):
        normal=cross(diff(pts[j],pts[i]),diff(pts[k],pts[i]))
        sides={v:dot(normal,diff(pts[v],pts[i])) for v in nodes}
        low,high=min(sides.values()),max(sides.values())
        if low==high or low<0<high:
            continue
        face=tuple(v for v in nodes if sides[v]==0)
        if face in facets:
            continue
        if high>0:
            normal=tuple(-x for x in normal)
        orbit=[min(face)]
        while len(orbit)<len(face):
            at=orbit[-1]
            next_point=next(v for v in face if v!=at and all(
                dot(normal,cross(diff(pts[v],pts[at]),diff(pts[w],pts[at])))>=0 for w in face))
            orbit.append(next_point)
        assert len(set(orbit))==len(face)
        facets[face]=tuple(orbit)
    return list(facets.values())

def lifted(face,negative):
    val=acb(1)
    for i,j in zip(face,face[1:]+face[:1]):
        pair=tuple(sorted((i,j)))
        ds,z=edges[pair]
        if ds==200 and pair not in negative:
            z=-z
        val*=z if i<j else z.conjugate()
    assert not val.imag.contains(0) or val.real>0
    return val.arg()

def degree(faces,negative):
    angles=[lifted(f,negative) for f in faces]
    # Each edge occurs with opposite orientation, so the sum is an integer multiple of 2pi.
    incidences={}
    for face in faces:
        for i,j in zip(face,face[1:]+face[:1]):
            incidences[i,j]=incidences.get((i,j),0)+1
    assert all(k==1 and incidences.get((j,i))==1 for (i,j),k in incidences.items())
    n=(sum(angles)/(2*arb.pi())).unique_fmpz()
    assert n is not None
    return int(n),min((arb.pi()-abs(a)).lower() for a in angles)

small,big=boundary(inner),boundary(wall)
signs=[pair for pair,(d,_) in edges.items() if d==200 and set(pair)<=set(wall)]
assert len(small)==4 and sorted(map(len,big))==[3]*4+[4]*3 and len(signs)==3
for mask in product((0,1),repeat=3):
    active={p for p,v in zip(signs,mask) if v}
    assert degree(small,active)[0]==len(active)//2
    assert degree(big,active)[0]==len(active)
active=set(signs)
margin=degree(big,active)[1]
quadrilaterals=[i for i,f in enumerate(big) if len(f)==4]
for choices in product((0,1),repeat=3):
    triangles=[]
    for i,f in enumerate(big):
        if len(f)==3:
            triangles.append(f)
            continue
        if choices[quadrilaterals.index(i)]:
            f=f[1:]+f[:1]
        triangles.extend(((f[0],f[1],f[2]),(f[0],f[2],f[3])))
    d,gap=degree(triangles,active)
    assert d==3
    margin=min(margin,gap)
print('wall vertices',wall,'inner vertices',inner,'negative edges',signs)
print('flux degree on inner and full wall hulls',degree(small,active)[0],degree(big,active)[0])
print('minimum principal branch margin',margin.str(15))
assert margin>arb('.0076')
