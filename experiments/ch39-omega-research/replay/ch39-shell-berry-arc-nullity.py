"""Flat components of the shortest-arc graphs, with exact spin-half Berry holonomy."""
from flint import fmpz_poly
from itertools import combinations, product
from math import isqrt


def berry_components(radius_squared):
    # A complex number is a pair of real polynomials evaluated at sqrt(R).
    r = radius_squared
    one = (fmpz_poly([1]),fmpz_poly([]))
    modulus = fmpz_poly([-r,0,1])
    pts = [p for p in product(range(isqrt(r)+1),repeat=3) if sum(x*x for x in p)==r]

    def mul(z,w):
        a,b = z
        c,d = w
        return ((a*c-b*d)%modulus,(a*d+b*c)%modulus)

    def conj(z):
        return (z[0], -z[1])

    def positive(p):
        a,b = int(p[0]),int(p[1])
        if not b:
            return a>0
        if not a or (a>0) == (b>0):
            return b>0
        return (a*a>r*b*b) == (a>0)

    edges = {}
    for i,j in combinations(range(len(pts)),2):
        x,y,z = pts[i]
        u,v,w = pts[j]
        cost = r-x*u-y*v-z*w
        ov = (fmpz_poly([r+z*w+x*u+y*v,z+w]),fmpz_poly([x*v-y*u]))
        edges.setdefault(cost,[]).append((i,j,ov))
    adjacency = [[] for p in pts]
    for cost,added in sorted(edges.items()):
        for i,j,ov in added:
            adjacency[i].append((j,ov))
            adjacency[j].append((i,conj(ov)))
        paths = {}
        components = flat = 0
        for vertex in range(len(pts)):
            if vertex in paths:
                continue
            components += 1
            paths[vertex] = one
            queue = [vertex]
            trivial = True
            for i in queue:
                for j,ov in adjacency[i]:
                    along = mul(paths[i],ov)
                    if j not in paths:
                        paths[j] = along
                        queue.append(j)
                    elif trivial:
                        cycle = mul(along,conj(paths[j]))
                        trivial = cycle[1].is_zero() and positive(cycle[0])
            flat += int(trivial)
        yield cost,components,flat


if __name__ == "__main__":
    expected = {17:(0,3),26:(0,3),29:(0,3),867:(1,76)}
    for r,(limit,first_cutoff) in expected.items():
        rows = []
        for cutoff,ncomp,flat in berry_components(r):
            rows.append((cutoff,ncomp,flat))
            if flat <= limit:
                break
        assert rows[-1][0] == first_cutoff and rows[-1][2] == 0
        print("squared radius",r,"(cutoff, components, parallel sections)",rows)
