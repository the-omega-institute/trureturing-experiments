"""Surface Haar integrals and an exact three-channel crossing on a finite complex."""
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import json

PLANE = [(1, 2, 3), (1, 2, 4), (1, 3, 5), (1, 4, 6), (1, 5, 6),
         (2, 3, 6), (2, 4, 5), (2, 5, 6), (3, 4, 5), (3, 4, 6)]
DIMS = (6, 3, 3, 2, 1)
TARGET = F(4322773606174, 3962837583051)


def topology():
    edges = Counter(edge for face in PLANE for edge in combinations(face, 2))
    assert len(edges) == 15 and set(edges.values()) == {2}
    for vertex in range(1, 7):
        graph = defaultdict(set)
        for face in PLANE:
            if vertex in face:
                a, b = (v for v in face if v != vertex)
                graph[a].add(b)
                graph[b].add(a)
        assert len(graph) == 5 and {len(x) for x in graph.values()} == {2}
        reached = {next(iter(graph))}
        for _ in graph:
            reached |= {v for u in reached for v in graph[u]}
        assert reached == set(graph)
    assert 6-len(edges)+len(PLANE) == 1
    boundary = PLANE[0]
    for edge in combinations(boundary, 2):
        edges[edge] -= 1
    assert {e for e, count in edges.items() if count == 1} == set(combinations(boundary, 2))


def orientable_sheet(n):
    q, u, d, l, e = (F(1, dim**n) for dim in DIMS)
    return 3*(q+4*u+d+l+e)/5, 9*q+l, 6*q+3*(u+d)/2


def disk_slope(n):
    a, b, c = orientable_sheet(n)
    return (a-b)/(b-c)


def cap(z):
    # A three-face disk with product of its face strengths equal to z.
    return tuple(v*z/(1+F(14, 3)*z) for v in (F(71, 45), F(1, 3), F(5, 18)))


def crosscap(t, s):
    # Each of the nine triangles has its own three-face cap. Haar coefficient
    # after the crosscap gluing is the ninth power of its Fourier eigenvalue.
    h = t*s**3
    vac, weak, color, both = (1+3*h)**9, (5*h/27)**9, (h/16)**9, (h/432)**9
    den = vac+3*weak+8*color+24*both
    return F(0), (2*weak+16*both)/den, (3*color+9*both)/den


def match(t):
    zero, w, c = crosscap(t, t)
    assert not zero and w > 2*c > 0
    z = c/(F(13, 10)-F(14, 3)*c)
    caps = cap(z)
    external = [a+b for a, b in zip(caps, (zero, w, c), strict=True)]
    assert external[0] == external[2] and external[1] > TARGET*external[0]
    common = (external[1]-TARGET*external[0])/(TARGET-1)
    t0 = common/(4-30*common)
    tuned = z/t**2
    assert all(0 < value < F(1, 15) for value in (t, tuned, t0))
    total = [v+4*t0/(1+30*t0) for v in external]
    assert total[0] == total[2] and total[1] == TARGET*total[2]
    return tuned, t0, total[2]


def verify():
    topology()
    assert [disk_slope(n) for n in (3, 5, 7, 9)] == [F(112, 5), F(776, 25), F(19936, 215), F(1006712, 3025)]
    y, w, c = cap(F(1))
    assert (y-w-F(218, 115)*(w-c))*(1+F(14, 3)) == F(131, 115)
    t = F(1, 30)
    for label in ('1/30', 'first archived positive zero'):
        if label != '1/30':
            data = json.loads(Path(__file__).with_name('ch39-tick-separators.json').read_text())
            cell = data['prefix'][0][0]
            height = F(2*cell+1, 2**401)
            t = 1/(F(1, 4)+height**2)
        values = match(t)
        print(label, 'two fitted strengths, common color/hypercharge stiffness:',
              [format(float(v), '.8e') for v in values])


if __name__ == '__main__':
    verify()
