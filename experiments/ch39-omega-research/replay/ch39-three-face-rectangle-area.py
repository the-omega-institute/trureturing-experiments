"""Single-slab center twists, cross-edge incidences and rectangular Wilson bounds."""
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations, permutations, product


def slab_incidences():
    origin = (0,0,0,0)
    faces = set()
    for corner in product((-1,0),repeat=4):
        for perm in permutations(range(4)):
            vertex = list(corner)
            simplex = [tuple(vertex)]
            for j in perm:
                vertex[j] += 1
                simplex.append(tuple(vertex))
            for face in combinations(simplex,3):
                if origin in face and {v[0] for v in face}=={0,1}:
                    faces.add(face)
    incidence = defaultdict(set)
    for face in faces:
        cross = [edge for edge in combinations(face,2) if edge[0][0]!=edge[1][0]]
        assert len(cross)==2
        for edge,neighbor in (cross,cross[::-1]):
            if origin in edge:
                incidence[edge].add(neighbor)
    assert len(incidence)==8 and max(map(len,incidence.values()))==14


def center_strip():
    # Pin both horizontal time slices. Sum the three vertical Z5 links of a width-two strip.
    weights = (6,5,5,5,5)
    for h,k in product(range(5),repeat=2):
        endpoints = Counter()
        for a,b,c in product(range(5),repeat=3):
            endpoints[a,c] += weights[(a-b+h)%5]*weights[(b-c+k)%5]
        assert all(endpoints[a,c]==endpoints[(a+2)%5,(c+2)%5] for a,c in endpoints)
        assert all(endpoints[a,c] == 135+int((a-c+h+k)%5==0) for a,c in endpoints)
    slab_rate = F(245,2000)
    # |1-zeta^2| > 19/10, so the bound for one unit of width is below 0.147.
    assert 2*slab_rate/((1-slab_rate)*F(19,10)) < F(147,1000)
    transported = F(139,1000)
    assert slab_rate**2*F(9,7) < transported**2
    assert 2*transported/((1-transported)*F(19,10)) < F(171,1000)


if __name__ == "__main__":
    slab_incidences()
    center_strip()
    print("Single-layer incidence bound 14, center twist, and rectangular area factor < 0.147.")
