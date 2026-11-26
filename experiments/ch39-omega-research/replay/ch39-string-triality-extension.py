"""Order-three defect class, cyclic lift, and the minimum set over the shell."""
from collections import Counter
from itertools import permutations, product


def cocycle(a,b,c):
    return (a%3)*((b%3+c%3)//3)%3


def check_lift():
    def beta(a,b):
        return (a%3)*(b//3)%3
    for a,b,c in product(range(9),repeat=3):
        cob = beta(b,c)-beta((a+b)%9,c)+beta(a,(b+c)%9)-beta(a,b)
        assert cob%3 == cocycle(a,b,c)
    assert all(3*x%9 for x in (1,4,7))  # no lifted generator for a Z3 lens holonomy


def check_shell():
    pts = {p for p in product(range(30),repeat=3) if sum(x*x for x in p)==867}
    def rotate(p):
        return p[1],p[2],p[0]
    def dist(p,q):
        return sum((x-y)**2 for x,y in zip(p,q))
    assert len(pts)==19
    pole = (17,17,17)
    near = sorted(p for p in pts if dist(p,pole)==68)
    assert len(near)==3 and set(near)=={(13,13,23),(13,23,13),(23,13,13)}
    assert len({tuple(dist(p,q) for q in [pole]+near) for p in pts})==19
    assert all(sum(dist(p,q)==68 for q in pts)==3
               for p in pts if len(set(p))==1)
    fingerprints = Counter(tuple(sorted(dist(p,q) for q in pts)) for p in pts)
    assert sorted(fingerprints.values())==[1,3,3,6,6]
    assert fingerprints[tuple(sorted(dist(pole,q) for q in pts))]==1
    todo = set(pts)
    comps = []
    while todo:
        part = {todo.pop()}
        stack = list(part)
        while stack:
            p = stack.pop()
            others = {q for q in todo if dist(p,q)<=56}
            todo -= others
            part |= others
            stack.extend(others)
        comps.append(part)
    assert Counter(map(len,comps))=={2:3,4:3,1:1}
    for comp in comps:
        if len(comp)==4:
            seq = list(comp)
            autos = [perm for perm in permutations(seq) if
                     all((dist(p,q)<=56)==(dist(u,v)<=56)
                         for p,u in zip(seq,perm) for q,v in zip(seq,perm))]
            assert len(autos)==4
    chosen = [near[0],rotate(near[0]),rotate(rotate(near[0]))]
    cover = {(p,k) for p in pts for k in range(3 if p in chosen else 1)}
    def lift(pair):
        p,k = pair
        return rotate(p), (k+int(p==chosen[-1]))%3
    lengths = []
    for v in cover:
        w = lift(v)
        n = 1
        while w!=v and n<10:
            w = lift(w)
            n += 1
        assert w==v and lift(v)[0]==rotate(v[0])
        lengths.append(n)
    assert len(cover)==25 and Counter(lengths)=={1:1,3:15,9:9}
    print("Z9 two-cochain verified on 729 triples; shell lifts from 19 to 25 states.")


if __name__ == "__main__":
    check_lift()
    check_shell()
