"""Shortest synchronous charged transfers for the stationary three-arrival filter."""
from collections import Counter
from itertools import product
from pathlib import Path
from runpy import run_path

root = Path(__file__).parent
shell = run_path(str(root/'ch39-local-shell-return.py'))
sites = [(1,5,29),(1,29,5),(7,17,23),(7,23,17),(11,11,25),(13,13,23)]
per_site = [36,36,36,36,16,48]


def orbit(p):
    rotate = shell['rotate']
    return [shell['scaled'](q,3) for q in (p,rotate(p),rotate(rotate(p)))]


def routes():
    d2 = shell['squared']
    adjacency = {p: dict((q,w) for q,w in row if d2(p,q)>=108)
                 for p,row in shell['graph'](3,56).items()}
    pole = (51,51,51)
    adjacency[pole] = {p:0 for p in adjacency if 108<=d2(p,pole)<=504}
    for p in adjacency[pole]:
        adjacency[p][pole] = 0
    assert len(adjacency)==76
    rotate = shell['rotate']
    assert all(adjacency[rotate(p)][rotate(q)]==w for p,row in adjacency.items() for q,w in row.items())
    layers = [[{(orbit(s)[0],0): (orbit(s)[0],)} for s in sites]]
    for _ in range(5):
        step = []
        for prev in layers[-1]:
            following = {}
            for (p,total), path in prev.items():
                for q,w in adjacency[p].items():
                    following.setdefault((q,total+w),path+(q,))
            step.append(following)
        layers.append(step)
    for length, sources in [(1,[0,1]),(2,[0,1]),(3,[0,1]),(4,[0,1,3])]:
        targets = {j for i in sources for j,s in enumerate(sites)
                   if any((q,1) in layers[length][i] for q in orbit(s))}
        assert sum(per_site[j] for j in targets) < sum(per_site[i] for i in sources)
    offsets = {}
    for i,j in product(range(6),repeat=2):
        if i<2 and j<2:
            continue
        table = layers[5][i]
        shift = next(d for d,q in enumerate(orbit(sites[j])) if all((q,k) in table for k in (-1,0,1)))
        end = orbit(sites[j])[shift]
        offsets[i,j] = shift
        for k in (-1,0,1):
            trace = table[end,k]
            assert len(trace)==6 and trace[0]==orbit(sites[i])[0] and trace[-1]==end
            assert sum(adjacency[a][b] for a,b in zip(trace,trace[1:]))==k
    return offsets


def position(color, digit):
    if digit < 4:
        return 4 if color < 4 else 5
    return color//4


def verify():
    offsets = routes()
    assert Counter(position(h,a) for h,a in product(range(16),range(13)))==dict(enumerate(per_site))
    def advance(state):
        o,h,a,b,c = state
        i,j = position(h,a),position(h^8,b)
        return (o+offsets[i,j])%3, h^8, b, c, a
    counts = Counter()
    images = set()
    for state in product(range(3),range(16),range(13),range(13),range(13)):
        nxt = advance(state)
        images.add(nxt)
        later = advance(nxt)
        counts[tuple(int(item[2]<4) for item in (state,nxt,later))]+=1
    assert len(images)==624*169
    for word in product(range(2),repeat=3):
        n = sum(word)
        assert counts[word]==48*4**n*9**(3-n)
    print('5 local edges per arrival, 32 ordered orbit pairs, three charges per pair.')
    print('Bijection on',len(images),'endpoint states; word counts',sorted(counts.items()))


if __name__=='__main__':
    verify()
