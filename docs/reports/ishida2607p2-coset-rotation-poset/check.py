#!/usr/bin/env python3
"""Independent exact checks for Ishida 2607.10306, Problem 2.
No third-party code or numeric approximations. Run with Python >= 3.10.
"""
from itertools import product, permutations
from random import Random
from collections import defaultdict
from hashlib import sha256
import json


def table(elements, op):
    ix = {g: i for i, g in enumerate(elements)}
    return [[ix[op(g, h)] for h in elements] for g in elements]


def cyclic(n):
    return table(list(range(n)), lambda x, y: (x+y) % n)


def abelian(*ns):
    es = list(product(*(range(n) for n in ns)))
    return table(es, lambda x, y: tuple((a+b) % n for a, b, n in zip(x,y,ns)))


def dihedral(k):
    es = list(product(range(k), range(2)))
    return table(es, lambda x,y: ((x[0]+(-1 if x[1] else 1)*y[0]) % k, (x[1]+y[1]) % 2))


def quaternion():
    # Elements (-1)^s * (1,i,j,k)[a], with i*j=k, j*k=i, k*i=j.
    es = list(product(range(2),range(4)))
    def op(x,y):
        s,a=x; t,b=y
        if a==0: return ((s+t)%2,b)
        if b==0: return ((s+t)%2,a)
        if a==b: return ((s+t+1)%2,0)
        return ((s+t+(0 if (a,b) in ((1,2),(2,3),(3,1)) else 1))%2,6-a-b)
    return table(es,op)


def validate(T):
    n=len(T)
    assert all(sorted(row)==list(range(n)) for row in T)
    assert all(T[0][a]==T[a][0]==a for a in range(n))
    assert all(T[T[a][b]][c]==T[a][T[b][c]] for a,b,c in product(range(n),repeat=3))
    return [next(b for b in range(n) if T[a][b]==T[b][a]==0) for a in range(n)]


def preferences(T,A):
    n=len(T); inv=validate(T)
    men=[[T[g][a] for a in A] for g in range(n)]
    women=[[T[h][inv[a]] for a in reversed(A)] for h in range(n)]
    rm=[[0]*n for _ in range(n)]; rw=[[0]*n for _ in range(n)]
    for g in range(n):
        for i,h in enumerate(men[g]): rm[g][h]=i
    for h in range(n):
        for j,g in enumerate(women[h]): rw[h][g]=j
    assert all(rm[g][h]+rw[h][g]==n-1 for g,h in product(range(n),repeat=2))
    return inv,men,rm,rw


def stable_matchings(T,A,men,rm,rw):
    n=len(T); out=[]
    for mu in permutations(range(n)):
        owner=[0]*n
        for g,h in enumerate(mu): owner[h]=g
        blocked=False
        for g in range(n):
            for h in men[g][:rm[g][mu[g]]]:
                if rw[h][g]<rw[h][owner[h]]:
                    blocked=True;break
            if blocked:break
        if not blocked:out.append(mu)
    return set(out)


def coset_poset(T,A,inv):
    n=len(T); nodes=[]
    for b in range(1,n):
        q=T[A[b-1]][inv[A[b]]]; H={0};x=q
        while x not in H: H.add(x);x=T[x][q]
        cs={frozenset(T[g][h] for h in H) for g in range(n)}
        for C in sorted(cs,key=lambda s:tuple(sorted(s))):nodes.append((b,C))
    pred=[]
    for j,(b,D) in enumerate(nodes):
        m=0
        for i,(a,C) in enumerate(nodes[:j]):
            if a+1==b and C&D:m|=(1<<i)|pred[i]
        pred.append(m)
    return nodes,pred


def ideal_matchings(T,A,nodes,pred):
    n=len(T); maps={}
    def rec(j,mask):
        if j==len(nodes):
            f=[sum(bool(mask>>i&1) and g in C for i,(b,C) in enumerate(nodes)) for g in range(n)]
            mu=tuple(T[g][A[f[g]]] for g in range(n))
            assert len(set(mu))==n
            assert mu not in maps
            maps[mu]=mask;return
        rec(j+1,mask)
        if pred[j]&mask==pred[j]:rec(j+1,mask|(1<<j))
    rec(0,0)
    return maps


def exposed_rotations(mu,men,rm,rw):
    n=len(mu); owner=[0]*n
    for g,h in enumerate(mu):owner[h]=g
    nxt={}
    for g in range(n):
        for h in men[g][rm[g][mu[g]]+1:]:
            if rw[h][g]<rw[h][owner[h]]:nxt[g]=owner[h];break
    seen=set();out={}
    for root in range(n):
        if root in seen:continue
        path=[];where={};x=root
        while x in nxt and x not in seen and x not in where:
            where[x]=len(path);path.append(x);x=nxt[x]
        if x in where:
            cyc=path[where[x]:];assert len(cyc)>=2
            rho=frozenset((g,mu[g],mu[nxt[g]]) for g in cyc)
            nu=list(mu)
            for g,old,new in rho:nu[g]=new
            out[rho]=tuple(nu)
        seen.update(path)
    return out


def check(T,A):
    inv,men,rm,rw=preferences(T,A);n=len(T)
    stable=stable_matchings(T,A,men,rm,rw)
    nodes,pred=coset_poset(T,A,inv);ideals=ideal_matchings(T,A,nodes,pred)
    assert stable==set(ideals), ('stable-vs-ideals',T,A)
    # Source Lemma 5.5 and Theorem 5.8, independently checked on brute-force stable states.
    for mu in stable:
        f=[rm[g][mu[g]] for g in range(n)];b=max(f)
        if b:
            q=T[A[b-1]][inv[A[b]]]
            assert all(f[T[g][q]]==b for g in range(n) if f[g]==b)
    lower=n+sum((1<<sum(a==b for a,C in nodes))-2 for b in range(1,n))
    assert len(stable)>=lower
    # The rotation computation above never receives Q or an ideal.
    edges={mu:exposed_rotations(mu,men,rm,rw) for mu in stable}
    rotations=set().union(*(set(es) for es in edges.values()))
    canonical={frozenset((g,T[g][A[b-1]],T[g][A[b]]) for g in C):i for i,(b,C) in enumerate(nodes)}
    assert len(canonical)==len(nodes) and rotations==set(canonical)
    start=tuple(men[g][0] for g in range(n)); eliminated={start:frozenset()}
    for mu in sorted(stable,key=lambda m:sum(rm[g][m[g]] for g in range(n))):
        assert mu in eliminated
        for rho,nu in edges[mu].items():
            assert nu in stable and rho not in eliminated[mu]
            new=eliminated[mu]|{rho}
            if nu in eliminated:assert eliminated[nu]==new
            else:eliminated[nu]=new
    assert len(eliminated)==len(stable)
    for rho in rotations:
        for sigma in rotations:
            # rho precedes sigma iff rho was eliminated before every exposure/elimination of sigma.
            precedes=rho==sigma or all(rho in eliminated[mu] for mu in stable if sigma in edges[mu])
            i,j=canonical[rho],canonical[sigma]
            assert precedes==(i==j or bool(pred[j]>>i&1))
    for mu in stable:
        assert {canonical[r] for r in eliminated[mu]}=={i for i in range(len(nodes)) if ideals[mu]>>i&1}
        for nu in stable:
            dominates=all(rm[g][mu[g]]<=rm[g][nu[g]] for g in range(n))
            assert dominates==(ideals[mu]&ideals[nu]==ideals[mu])
    return len(stable),len(nodes),lower


def main():
    groups=[('C1',cyclic(1)),('C2',cyclic(2)),('C3',cyclic(3)),('C4',cyclic(4)),
            ('C2xC2',abelian(2,2)),('C5',cyclic(5)),('C6',cyclic(6)),('S3',dihedral(3)),
            ('C7',cyclic(7)),('C8',cyclic(8)),('C4xC2',abelian(4,2)),('C2xC2xC2',abelian(2,2,2)),
            ('D8',dihedral(4)),('Q8',quaternion())]
    rng=Random(260710306);total=0;records=[]
    for name,T in groups:
        n=len(T)
        if n<=5: orders=list(permutations(range(n)))
        else:
            pool={tuple(range(n)),tuple(reversed(range(n)))}
            if name=='C6':pool.add((0,2,4,1,3,5))
            while len(pool)<40:
                a=list(range(n));rng.shuffle(a);pool.add(tuple(a))
            orders=sorted(pool)
        counts=[]
        for A in orders:
            z,r,lower=check(T,A);counts.append(z);total+=1
            if name=='C2xC2' and A==(0,2,1,3):
                assert (z,r)==(10,6);print('FIDELITY Klein source ordering ((0,0),(1,0),(0,1),(1,1)): stable=10 rotations=6')
            if name=='C4' and A==(0,2,1,3):
                assert z==8;print('FIDELITY C4 ordering (0,2,1,3): stable=8')
            if name=='C6' and A==(0,2,4,1,3,5):
                assert (z,r)==(24,11);print('SANITY C6 ordering (0,2,4,1,3,5): stable=24 rotations=11')
            records.append((name,A,z,r,lower))
        print(f'{name}: order={n} orderings={len(orders)} stable_min={min(counts)} stable_max={max(counts)}')
    digest=sha256(json.dumps(records,separators=(',',':')).encode()).hexdigest()
    print(f'Exact cases={total}; groups={len(groups)}; all orderings through order 5; 40 orderings per larger group')
    print('Per-case integer results sha256='+digest)
    print('Checks: source rank-sum, top-coset lemma, counting lower bound, stable bijection, male order, independent exposed rotations and precedence')
    print('ALL_CHECKS_OK')


if __name__=='__main__':main()
