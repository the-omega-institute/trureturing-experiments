"""度量三角驻点的剩余低维 Hessian 与零通量度量回路的刚化。"""
from collections import defaultdict
from itertools import combinations
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat, arb_poly

phase = run_path(str(Path(__file__).with_name("ch39-matter-number-three-spokes.py")))
soft = phase["soft"]
pole = soft["pole"]
orbits = soft["star"]["charged"]["orbits"]
active, mirror = map(frozenset, orbits[:2])


def free_curvature(matrices, weight, c, moved):
    sine = (1-c*c).sqrt()
    result = arb(0)
    for q, matrix in matrices.items():
        values = []
        for t in (acb(1), acb(-1), acb(0,1), acb(0,-1)):
            m = acb_mat(matrix)
            for j in moved:
                m[pole,j], m[j,pole] = -t, -t.conjugate()
            a,b,d = phase["spoke_polynomials"](m,pole,orbits[0],weight,acb(1))
            root = acb(c,sine)
            values.append(a+b*root.real+d*root.imag)
        start, opposite, left, right = values
        # A common shift on any chosen spokes has Fourier degree one.
        cosine, linear = (start-opposite)*arb(".5"), (left-right)*arb(".5")
        for p in (cosine, linear):
            assert all(p[k].contains(0) for k in (17,18,19))
        cosine, linear = (arb_poly(p.coeffs()[:17]) for p in (cosine,linear))
        numerator = cosine*start+linear**2
        result += {0:33,5:18,10:3}[abs(q)]*phase["integrate_force"](numerator,start**2,"2e-5")
    return result/(8*arb.pi()**2)


def classify_faces():
    geo = soft["geo"]
    lengths = {(i,j):d for i,j,d,_ in geo["edges"]}
    vertices = set(range(len(geo["sites"])))-{pole}
    shapes = defaultdict(list)
    for i,j in combinations(sorted(vertices),2):
        for path in [(i,j), *((i,k,j) for k in sorted(vertices-{i,j}))]:
            walk = (pole,*path,pole)
            ds = tuple(lengths.get(tuple(sorted(pair))) for pair in zip(walk,walk[1:]))
            if None not in ds:
                shapes[min(ds,ds[::-1])].append((i,j))
    components = [{v} for v in vertices]
    for pairs in shapes.values():
        if any((i in active) != (j in active) for i,j in pairs):
            continue
        for i,j in pairs:
            incident = [part for part in components if i in part or j in part]
            components = [part for part in components if part not in incident]+[set.union(*incident)]
    assert set(map(frozenset,components)) == {active,mirror,frozenset(vertices-active-mirror)}
    return shapes


def verify():
    matrices,weight,c,weights,signatures = phase["verify"]()
    shapes = classify_faces()
    coeff = dict(zip(signatures,weights))
    def second(moved):
        value = free_curvature(matrices, weight, c, moved)
        for signature,pairs in shapes.items():
            if len(signature) != 3:
                continue
            k = coeff.get(tuple(sorted(signature)),0)
            for i,j in pairs:
                if (i in moved) != (j in moved):
                    value += k*(c if (i in active) != (j in active) else 1)
        return value
    h1,h2,h12 = (second(moved) for moved in (active,mirror,active|mirror))
    cross = (h12-h1-h2)/2
    determinant = h1*h2-cross**2
    print("Two orbit angles modulo bulk shift: diagonal, cross, determinant",
          *(v.str(8) for v in (h1,h2,cross,determinant)))
    assert h1 > arb(".017") and determinant > arb(".0024")

    sites = sorted(j for group in orbits for j in group)
    rotation = soft["star"]["charged"]["rotation"]
    cache = {}
    def direction(vertices):
        moved = frozenset(vertices)
        twice = frozenset(rotation[i] for i in moved)
        thrice = frozenset(rotation[i] for i in twice)
        key = min(tuple(sorted(v)) for v in (moved,twice,thrice))
        if key not in cache:
            cache[key] = second(frozenset(key))
        return cache[key]
    # Work in the gauge where the last spoke angle is held fixed at zero.
    basis = sites[:-1]
    size = len(basis)
    full = arb_mat(size,size)
    safe = arb_mat(size,size)
    for i,j in enumerate(basis):
        full[i,i] = direction([j])
    for i,j in combinations(range(size),2):
        v = (direction([basis[i],basis[j]])-full[i,i]-full[j,j])/2
        full[i,j] = full[j,i] = v
    for pairs in shapes.values():
        if any((i in active) != (j in active) for i,j in pairs):
            continue
        for i,j in pairs:
            for site in (i,j):
                if site in basis:
                    k = basis.index(site)
                    safe[k,k] += 1
            if i in basis and j in basis:
                k,l = basis.index(i),basis.index(j)
                safe[k,l] -= 1
                safe[l,k] -= 1

    def positive(weight):
        matrix = full+weight*safe
        mid = matrix.mid()
        # Weyl's bound from the absolute row sums of the symmetric radius.
        error = max(sum(matrix[i,j].rad() for j in range(size)).upper() for i in range(size))
        for k in range(size):
            mid[k,k] -= error
        return all(arb_mat([[mid[i,j] for j in range(k)] for i in range(k)]).det() > 0
                   for k in range(1,size+1))
    witness = arb_mat([[v] for v in (36,27,27,36,-61,-44,7,6,-3,5,5,5,5,8,8,-3,-3)])
    lower,upper = arb(".023"),arb(".0245")
    direction_value = (witness.transpose()*(full+lower*safe)*witness)[0,0]
    assert direction_value < 0 and positive(upper)
    print("Critical extra equal face weight bracket",lower.str(5),upper.str(5),
          "lower-end negative quadratic form upper bound",direction_value.upper().str(7))


if __name__ == "__main__":
    verify()
