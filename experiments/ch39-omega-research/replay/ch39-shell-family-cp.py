"""Test the antiunitary symmetry of family textures made from one or two arithmetic shells."""
from pathlib import Path
from itertools import combinations
from flint import acb,acb_mat,arb
import runpy

book=runpy.run_path(str(Path(__file__).with_name("ch39-book-shell-clock-embedding.py")))
cross,dot=book["cross"],book["dot"]
def unit(u):
    v=arb(dot(u,u)).sqrt()
    return [x/v for x in u]
def J(n):
    x,y,z=n
    return acb_mat([[0,-z,y],[z,0,-x],[-y,x,0]])*acb(0,1)
def family(sites,offset,c2num,c2den):
    normals=[unit(cross(sites[0],sites[j])) for j in (7,8)]
    n,m=normals
    assert (dot(n,m)**2-arb(c2num)/c2den).contains(0)
    amplitudes=[(2*book["h"]*f).sin() for f in book["f"][offset:offset+2]]
    A,B=[J(u)*a for u,a in zip(normals,amplitudes)]
    return A,B
a,b=family(book["v"],0,81,98)
c,d=family(book["image"],2,8741547,12563572)
odd=(a*b*c).trace().imag
assert not odd.contains(0)
Y1,Y2=a+b*b,c+d*d
def mass(Y):
    U=Y*Y
    vals=sorted((x.real for x in U.eig()),key=lambda v:float(v.mid()))
    assert vals[0]>0 and all(x<y for x,y in zip(vals,vals[1:]))
    return U,vals
H,x=mass(Y1)
K,y=mass(Y2)
def vandermonde(vals):
    p=arb(1)
    for u,v in combinations(vals,2):
        p*=v-u
    return p
j=abs((H*K-K*H).det().imag)/(2*vandermonde(x)*vandermonde(y))
assert j>0 and j<1/(6*arb(3).sqrt())
I=acb_mat([[int(i==k) for k in range(3)] for i in range(3)])
def projectors(H,values):
    ps=[]
    for i,root in enumerate(values):
        p=I
        for k,z in enumerate(values):
            if i!=k:
                p=(p*(H-I*z))/(root-z)
        ps.append(p)
    return ps
px,py=projectors(H,x),projectors(K,y)
mix=[(px[i]*py[k]).trace().real for i,k in ((0,1),(0,2),(1,2))]
assert all(0<v<1 for v in mix)
print("odd trace of three adjoints",odd.str(13))
print("same-prescription cross-shell flavor determinant |J|",j.str(13))
print("|Vus|, |Vub|, |Vcb|",*[v.sqrt().str(12) for v in mix])
