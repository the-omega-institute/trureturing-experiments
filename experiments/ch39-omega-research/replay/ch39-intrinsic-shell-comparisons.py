"""Build the 624 internal channels from signed chord comparisons on the shell."""
from collections import Counter
from fractions import Fraction as F
from itertools import product


def squared(p, q):
    return sum((x-y)**2 for x,y in zip(p,q))


def rotate(p):
    return p[1], p[2], p[0]


def verify():
    shell = {p for p in product(range(30),repeat=3) if sum(c*c for c in p)==867}
    profiles = {p: Counter(squared(p,q) for q in shell if q != p) for p in shell}
    poles = [p for p, dist in profiles.items() if len(dist)==3]
    assert poles == [(17,17,17)]
    pole = poles[0]
    near2 = min(profiles[pole])
    probes = [p for p in shell if squared(p,pole)==near2]
    assert len(probes)==3 and near2==68
    first = min(probes)
    chosen = [first,rotate(first),rotate(rotate(first))]
    assert set(chosen)==set(probes)
    def lift(state):
        p, j = state
        return rotate(p), (j+int(p==chosen[2]))%3
    readings = {}
    for p in shell-{pole}:
        values = {F(squared(p,pole)+near2-squared(p,b),2) for b in probes}
        direct = {sum((x-o)*(y-o) for x,y,o in zip(p,b,pole)) for b in probes}
        assert values == direct
        signed = values | {-v for v in values}
        assert len(signed)==2*len(values)
        readings[p] = signed
    states = [(p,k) for p in readings for k in range(3 if p in probes else 1)]
    weights = {s:len(readings[s[0]]) for s in states}
    assert Counter(weights.values()) == {4:12,6:12}
    channels = {(s,a,b) for s in states for a,b in product(readings[s[0]], repeat=2)}
    assert len(channels)==624 and sum(w*w for w in weights.values())==624
    assert sum(len(readings[p])**2 for p in readings)==528
    assert sum(weights.values())==120 and sum(w**3 for w in weights.values())==3360
    # Isometric embedding and compressed reflection, evaluated for all fiber pairs.
    fibers = [1,*[w*w for w in weights.values()]]
    roots = [1,*weights.values()]
    for i, ni in enumerate(fibers):
        for j, nj in enumerate(fibers):
            entry = F(2*ni*nj,625)-(ni if i==j else 0)
            assert entry/(roots[i]*roots[j]) == F(2*roots[i]*roots[j],625)-int(i==j)
    todo = set(channels)
    cycles = Counter()
    while todo:
        s0,a,b = next(iter(todo))
        s = s0
        length = 0
        while (s,a,b) in todo:
            todo.remove((s,a,b))
            length += 1
            s = lift(s)
        assert s==s0
        cycles[length] += 1
    assert cycles == {3:160,9:16}
    assert F(3*cycles[3],624)==F(10,13)
    print("One-register incidence cycles",dict(cycles),"third-step return",F(10,13))


if __name__ == "__main__":
    verify()
