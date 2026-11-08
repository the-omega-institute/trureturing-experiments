#!/usr/bin/env python3
"""Original exact audit, Python standard library only; no floating arithmetic.
Knill arXiv:2605.30799v1, Sections 1.7 and 3.1. Finite evidence, no Lean proof.
"""
import collections
import itertools as it
import json


def mul(g, h):
    return tuple(g[h[i]] for i in range(len(g)))


def cyc(n, c):
    p = list(range(n))
    for a, b in zip(c, c[1:]+c[:1]):
        p[a] = b
    return tuple(p)


def parity(p):
    return sum(p[i] > p[j] for i in range(len(p)) for j in range(i+1,len(p))) % 2


def smul(g, h, k):
    p,a = g
    q,b = h
    return mul(p,q), tuple((b[i]+a[q[i]]) % k for i in range(len(p)))


def cmul(g,h):
    return smul(g[0],h[0],3), smul(g[1],h[1],2)


def legal(g):
    (cp,co),(ep,eo) = g
    return (sorted(cp) == list(range(8)) and sorted(ep) == list(range(12))
            and sum(co)%3 == sum(eo)%2 == 0 and parity(cp) == parity(ep))


def graph():
    a = [set() for _ in range(16)]
    for i in range(8):
        for x,y in [(i,(i+1)%8),(i,8+i),(8+i,8+(i+3)%8)]:
            a[x].add(y)
            a[y].add(x)
    assert all(len(x)==3 for x in a) and sum(map(len,a))//2==24
    return a


def distances(a):
    ds = []
    for v in range(16):
        d = [None]*16
        d[v] = 0
        q = collections.deque([v])
        while q:
            x = q.popleft()
            for y in a[x]:
                if d[y] is None:
                    d[y] = d[x]+1
                    q.append(y)
        assert None not in d
        ds.append(d)
    return ds


def automorphisms(d):
    """Exhaustive bijection search; prune only distance-incompatible images.
    Every automorphism preserves distances. Every leaf preserves adjacency=d=1.
    At each node ALL compatible unused images are visited, hence completeness.
    """
    used, free, answer, nodes = {}, set(range(16)), [], 0
    def visit():
        nonlocal nodes
        nodes += 1
        if len(used)==16:
            answer.append(tuple(used[i] for i in range(16)))
            return
        choices = [(len(w),v,w) for v in range(16) if v not in used
                   for w in [[x for x in sorted(free)
                              if all(d[v][i]==d[x][j] for i,j in used.items())]]]
        _,v,ws = min(choices)
        for w in ws:
            used[v] = w
            free.remove(w)
            visit()
            free.add(w)
            del used[v]
    visit()
    return answer,nodes


PAIRS = [(r,r+4) for r in (0,1,2,3,8,9,10,11)]
COORD = {v:(i,t) for i,p in enumerate(PAIRS) for t,v in enumerate(p)}


def signed(g):
    p,a = zip(*(COORD[g[x]] for x,_ in PAIRS))
    assert sorted(p)==list(range(8))
    assert all(g[v]==PAIRS[p[i]][(t+a[i])%2]
               for i,pair in enumerate(PAIRS) for t,v in enumerate(pair))
    return tuple(p),tuple(a)


def embed(s):
    p,a = s
    m = len(p)
    assert m<=11
    ep = tuple(p)+tuple(range(m,12))
    eo = tuple(a)+(sum(a)%2,)+(0,)*(11-m)
    cp = cyc(8,[0,1]) if parity(p) else tuple(range(8))
    return (cp,(0,)*8),(ep,eo)


def closure(gens,n):
    one = tuple(range(n))
    seen,todo = {one},[one]
    while todo:
        x = todo.pop()
        for g in gens:
            y = mul(g,x)
            if y not in seen:
                seen.add(y)
                todo.append(y)
    return seen


CC = ([0,1,2,3],[0,3,7,4],[0,4,5,1],[4,7,6,5],[1,5,6,2],[2,6,7,3])
EC = ([0,1,2,3],[0,11,4,8],[1,8,5,9],[4,7,6,5],[2,9,6,10],[3,10,7,11])
TW = ({},{0:2,3:1,4:1,7:2},{0:1,1:2,4:2,5:1},{},{1:1,2:2,5:2,6:1},{2:1,3:2,6:2,7:1})
FL = ((),(),(1,5,8,9),(),(),(3,7,10,11))


def face(f):
    cp,ep = cyc(8,CC[f]),cyc(12,EC[f])
    return ((cp,tuple(TW[f].get(cp[i],0) for i in range(8))),
            (ep,tuple(int(ep[i] in FL[f]) for i in range(12))))


def geometry():
    """Derive face turns independently by exact integer 3D sticker rotations."""
    ns = dict(U=(0,1,0),R=(1,0,0),F=(0,0,1),D=(0,-1,0),L=(-1,0,0),B=(0,0,-1))
    cs = ['URF','UFL','ULB','UBR','DFR','DLF','DBL','DRB']
    es = ['UR','UF','UL','UB','DR','DF','DL','DB','FR','FL','BL','BR']
    def pos(label):
        return tuple(sum(ns[c][i] for c in label) for i in range(3))
    def rot(x,n):
        a,b,c = n
        u,v,w = x
        d = a*u+b*v+c*w
        return a*d-(b*w-c*v),b*d-(c*u-a*w),c*d-(a*v-b*u)
    ans = []
    for f in 'URFDLB':
        n,parts = ns[f],[]
        for ls,k in [(cs,3),(es,2)]:
            locs = list(map(pos,ls))
            p,a = [],[]
            for i,label in enumerate(ls):
                active = sum(x*y for x,y in zip(locs[i],n))==1
                j = locs.index(rot(locs[i],n)) if active else i
                ts = [ls[j].index(next(c for c in ns
                                      if ns[c]==(rot(ns[s],n) if active else ns[s])))
                      for s in label]
                assert all(ts[t]==(t+ts[0])%k for t in range(k))
                p.append(j)
                a.append(ts[0])
            parts.append((tuple(p),tuple(a)))
        ans.append(tuple(parts))
    return ans


def main():
    a = graph()
    d = distances(a)
    layers = [[v for v in range(16) if d[0][v]==k] for k in range(5)]
    assert list(map(len,layers))==[1,3,6,5,1]
    z = tuple((v//8)*8+(v%8+4)%8 for v in range(16))
    assert all([w for w in range(16) if d[v][w]==4]==[z[v]] for v in range(16))
    r = tuple((v//8)*8+(v%8+1)%8 for v in range(16))
    s = tuple((1-v//8)*8+(3*(v%8))%8 for v in range(16))
    for g in [r,s]:
        assert all((j in a[i])==(g[j] in a[g[i]]) for i in range(16) for j in range(16))
    assert {g[0] for g in closure([r,s],16)}==set(range(16))
    aut,nodes = automorphisms(d)
    aset = set(aut)
    assert len(aut)==len(aset)==96
    signed_images = {g:signed(g) for g in aut}
    images = {g:embed(signed_images[g]) for g in aut}
    assert len(set(images.values()))==96 and all(map(legal,images.values()))
    for g,h in it.product(aut,repeat=2):
        gh = mul(g,h)
        assert gh in aset and mul(g,z)==mul(z,g)
        assert signed_images[gh]==smul(signed_images[g],signed_images[h],2)
        assert images[gh]==cmul(images[g],images[h])
    gens,sub = [],{tuple(range(16))}
    for g in aut:
        if g not in sub:
            gens.append(g)
            sub = closure(gens,16)
    assert sub==aset and all(legal(images[g]) for g in gens)
    faces = [face(f) for f in range(6)]
    assert faces==geometry()
    one = ((tuple(range(8)),(0,)*8),(tuple(range(12)),(0,)*12))
    for f in faces:
        assert legal(f) and cmul(cmul(f,f),cmul(f,f))==one and cmul(f,f)!=one
    count = 0
    for m in range(4):
        w = [(p,e) for p in it.permutations(range(m)) for e in it.product(range(2),repeat=m)]
        assert len({embed(x) for x in w})==len(w) and all(legal(embed(x)) for x in w)
        for x,y in it.product(w,repeat=2):
            assert embed(smul(x,y,2))==cmul(embed(x),embed(y))
            count += 1
    print(json.dumps(dict(vertices=16,edges=24,distance_layers_u0=layers,
        automorphisms=96,backtracking_nodes=nodes,distinct_legal_images=96,
        multiplication_pairs=9216,automorphism_generators=gens,generator_images_legal=True,
        face_turns_checked_against_integer_geometry=6,universal_map_small_cases_m=list(range(4)),
        small_case_pairs=count,
        reachability_boundary='Constraints checked; sufficiency is classical, not proved here.'),indent=2))
    print('ALL_CHECKS_OK')


if __name__=='__main__':
    main()
