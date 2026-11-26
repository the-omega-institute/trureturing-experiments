"""Tied-shell flag complexes, their mirror symmetry, and the minimum finite holonomy for degree three."""
from itertools import combinations, permutations,product
from pathlib import Path
from runpy import run_path
from collections import Counter,deque
from flint import arb

m=run_path(str(Path(__file__).with_name("ch39-shell-third-shell.py")))
dot,cross=m["dot"],m["cross"]
shells=m["fattest"](6)


def rotate(v):
    return v[1:]+v[:1]


def collapse(r):
    vs=sorted(shells[r])
    cut=m["shells"][r][0][-1]
    assert all(len(set(p))==3 for p in vs)
    def D(i,j):
        return r-dot(vs[i],vs[j])
    cells=[sum(all(D(i,j)<=cut for i,j in combinations(c,2)) for c in combinations(range(12),k))
           for k in range(1,8)]
    S=set(range(12))
    while True:
        nbr={i:{j for j in S if D(i,j)<=cut} for i in S}
        opts=[(i,j) for i in sorted(S) for j in sorted(S) if i!=j and nbr[i]<=nbr[j]]
        if not opts:
            break
        S.remove(opts[0][0])
    print("r, simplex counts, core vertices",r,cells,len(S))
    return [vs[j] for j in sorted(S)]


def parity(v):
    return (-1)**sum(x>y for x,y in combinations(v,2))


def main():
    core=collapse(26)
    assert len(collapse(29))==1 and set(core)==set(permutations((1,3,4)))
    axes={(v.index(3),parity(v)):v for v in core}
    assert all((26-dot(v,w)<=7)==(v.index(3)!=w.index(3)) for v,w in combinations(core,2))
    faces=[]
    for signs in product((-1,1),repeat=3):
        f=tuple(axes[k,s] for k,s in enumerate(signs))
        if signs[0]*signs[1]*signs[2]<0:
            f=f[::-1]
        faces.append(f)
    def sides(f):
        return tuple(zip(f,f[1:]+f[:1]))
    oriented=Counter(e for f in faces for e in sides(f))
    assert all(n==1 and oriented[y,x]==1 for (x,y),n in oriented.items())
    # All coordinate permutations act by proper abstract rotations on this core.
    cuts=[]
    phases=[]
    for f in faces:
        v,w,z=f
        num=dot(v,cross(w,z))
        den=arb(26).sqrt()*(26+dot(v,w)+dot(w,z)+dot(z,v))
        phi=(arb(num)/den).atan()
        assert not phi.contains(0) and abs(phi)<arb.pi()/12
        phases.append(phi)
        cuts.append(2 if phi<0 else 1)
        for p in permutations(range(3)):
            mapped=tuple(tuple(v[i] for i in p) for v in f)
            assert any(mapped==t[k:]+t[:k] for t in faces for k in range(3))
    assert sum(phases).contains(0) and cuts.count(2)==4
    # Construct a Z/4 edge cochain on a dual spanning tree, then symmetrize over C3.
    parent={0:None}
    queue=deque([0])
    while queue:
        i=queue.popleft()
        for x,y in sides(faces[i]):
            j=next(k for k,f in enumerate(faces) if (y,x) in sides(f))
            if j not in parent:
                parent[j]=(i,y,x)
                queue.append(j)
    assert len(parent)==8
    co={e:0 for e in oriented}
    rhs=cuts[:]
    for j in reversed(list(parent)[1:]):
        i,x,y=parent[j]
        co[x,y]=rhs[j]%4
        co[y,x]=-co[x,y]%4
        rhs[i]+=rhs[j]
    assert rhs[0]%4==0
    field={}
    for x,y in co:
        vals=(co[x,y],co[rotate(x),rotate(y)],co[rotate(rotate(x)),rotate(rotate(y))])
        field[x,y]=3*sum(vals)%4
    actual=[]
    for f,offset,phi in zip(faces,cuts,phases):
        assert sum(field[e] for e in sides(f))%4==offset
        principal=phi+arb.pi()*offset/2
        assert 0<principal<arb.pi()
        actual.append(principal)
    assert (sum(actual)/(2*arb.pi())).unique_fmpz()==3
    for e,z in field.items():
        assert field[e[1],e[0]]==(-z)%4 and z==field[rotate(e[0]),rotate(e[1])]
    print("bare oriented face phases",*[p.str(7) for p in phases])
    print("order-four principal degree",3,"cut margin",min(arb.pi()-p for p in actual).str(12),
          "nonzero edges in this witness",sum(z!=0 for z in field.values())//2)


    return core,faces,phases,field


if __name__=="__main__":
    main()
