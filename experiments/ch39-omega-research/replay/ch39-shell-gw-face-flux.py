"""Edge twist of the local octahedral kernel by the actual determinant-lift face holonomies."""
from fractions import Fraction
from itertools import combinations, permutations, product
from pathlib import Path
from runpy import run_path
from flint import arb,acb
import sympy as sp


folder=Path(__file__).resolve().parent
shell=run_path(str(folder/"ch39-shell-local-octahedral-chirality.py"))
linear=run_path(str(folder/"ch39-shell-berry-rh-matching.py"))
vertices,H,_=shell["make_operator"]()
zero=arb(0)
def ball(value):
    re,im=sp.expand_complex(value).as_real_imag()
    def cv(x):
        if x.is_Rational:
            return arb(int(x.p))/int(x.q)
        # At this normalization all site entries belong to Q(sqrt(2)).
        a,b=x.expand().coeff(sp.sqrt(2),0),x.expand().coeff(sp.sqrt(2),1)
        assert a.is_Rational and b.is_Rational and sp.expand(a+b*sp.sqrt(2)-x)==0
        return cv(a)+arb(2).sqrt()*cv(b)
    return acb(cv(re),cv(im))


def maps():
    points=list(permutations((1,3,4)))
    indexed={}
    for p in points:
        sign=(-1)**sum(x>y for x,y in combinations(p,2))
        pos=p.index(3)
        idx=vertices.index(tuple(sign if j==pos else 0 for j in range(3)))
        indexed[idx]=p
    faces=[]
    for signs in product((-1,1),repeat=3):
        verts=tuple(vertices.index(tuple(s if j==i else 0 for j in range(3))) for i,s in enumerate(signs))
        faces.append(verts if signs[0]*signs[1]*signs[2]>0 else verts[::-1])
    return indexed,faces


def edges(face):
    return list(zip(face,face[1:]+face[:1]))


def cochain(faces,turns):
    visited={0:None}
    queue=[0]
    for i in queue:
        for x,y in edges(faces[i]):
            j=next(j for j,f in enumerate(faces) if (y,x) in edges(f))
            if j not in visited:
                visited[j]=(i,y,x)
                queue.append(j)
    rhs=list(turns)
    field={}
    for j in reversed(queue[1:]):
        i,x,y=visited[j]
        field[x,y]=rhs[j]%1
        field[y,x]=(-rhs[j])%1
        rhs[i]+=rhs[j]
    assert rhs[0]%1==0
    assert all((sum(field.get(e,0) for e in edges(face))-v)%1==0 for face,v in zip(faces,turns))
    return field


def twisted(field):
    out=[]
    for i in range(12):
        row=[]
        for j in range(12):
            turn=field.get((j%6,i%6),Fraction(0))
            rot=acb(0,2*arb.pi()*turn.numerator/turn.denominator).exp()
            row.append(rot*ball(H[i,j]))
        out.append(row)
    return out


if __name__ == "__main__":
    indexed,faces=maps()
    offsets=[]
    for f in faces:
        ds=[26-sum(a*b for a,b in zip(indexed[i],indexed[j])) for i,j in edges(f)]
        offsets.append(Fraction(1,4) if ds==[7]*3 else Fraction(5,12))
    assert offsets.count(Fraction(1,4))==2
    for name,turns in (("balanced shell lifts",offsets),("trivial physical edge connection",[Fraction(0)]*8)):
        relative=[v-Fraction(1,8) for v in turns]
        matrix=twisted(cochain(faces,relative))
        count=linear["count_below"]
        neg=count(matrix,-arb(1)/100)
        assert neg==count(matrix,arb(1)/100)
        print(name,"negative inertia",neg,"net index",6-neg,"gap > 1/100")
