"""Certify a product box in the negative Regge tail and its normalized probability."""

from itertools import combinations, permutations, product
from math import factorial
from pathlib import Path
from runpy import run_path

from flint import arb, ctx, fmpq, fmpq_mat, fmpq_mpoly_ctx

EDGES = list(combinations(range(5), 2))
VARIABLES = fmpq_mpoly_ctx.get(("q", 10)).gens()
FACES = run_path(str(Path(__file__).with_name("ch39-regge-shape-support.py")))["FACES"]
GAPS = [(1973, 2000), (1917, 2000), (674, 694), (1953, 2000)]
BOX = [tuple(fmpq(x, 2000) for x in GAPS[j-i-1]) for i, j in EDGES]


def gram(edges):
    q = dict(zip(EDGES, edges))
    return [[(q[0, i] + q[0, j] - (q[min(i, j), max(i, j)] if i != j else 0)) / 2
             for j in range(1, 5)] for i in range(1, 5)]


def determinant_polynomial():
    matrix = gram(VARIABLES)
    det = 0
    for columns in permutations(range(4)):
        term = (-1)**sum(columns[i] > columns[j] for i, j in combinations(range(4), 2))
        for i, j in enumerate(columns):
            term *= matrix[i][j]
        det += term
    return det


def range_on_unit_box(polynomial):
    # Nonconstant even monomials range in [0,1], all others in [-1,1].
    low = high = fmpq(0)
    for powers, coefficient in polynomial.terms():
        if not any(powers):
            low += coefficient
            high += coefficient
        elif any(power % 2 for power in powers):
            low -= abs(coefficient)
            high += abs(coefficient)
        else:
            low += min(0, coefficient)
            high += max(0, coefficient)
    return arb(low).union(arb(high))


def determinant_extrema(det, derivatives):
    minimum = fmpq(1)
    for vertex in product(*BOX):
        matrix = gram(vertex)
        for k in range(1, 5):
            minor = fmpq_mat([row[:k] for row in matrix[:k]]).det()
            assert minor > 0
        minimum = min(minimum, minor)
    upper = [hi for lo, hi in BOX]
    assert all(p(*upper) > 0 for p in derivatives)
    maximum = det(*upper)
    assert minimum > fmpq(1, 16)
    return minimum, maximum


def face_enclosures(box, det, derivatives):
    centers = [(a+b)/2 for a,b in box]
    radii = [(b-a)/2 for a,b in box]
    coordinates = [c+r*x for c,r,x in zip(centers, radii, VARIABLES)]
    q = [arb(c, r) for c,r in zip(centers, radii)]
    d = range_on_unit_box(det.compose(*coordinates))
    de = [range_on_unit_box(p.compose(*coordinates)) for p in derivatives]
    grad = [arb(0) for _ in EDGES]
    action = arb(0)
    for i, j, a, b, c, valence in FACES:
        opposite = EDGES.index((i,j))
        triangle = [EDGES.index(edge) for edge in [(a,b), (a,c), (b,c)]]
        u, v, w = (coordinates[e] for e in triangle)
        area = range_on_unit_box((2*(u*v+u*w+v*w)-u*u-v*v-w*w)/16).sqrt()
        angle = (de[opposite]/(2*area*d.sqrt())).atan()+arb.pi()*(fmpq(2,valence)-fmpq(1,2))
        action += area * angle
        for s in range(3):
            e, f, g = (triangle[(s+t)%3] for t in range(3))
            grad[e] += angle*(q[f]+q[g]-q[e])/(16*area)
    return d, de, action, grad


def local_maximum(det, derivatives):
    box = BOX.copy()
    while any(lo < hi for lo, hi in box):
        centers = [((a+b)/2, (a+b)/2) for a,b in box]
        _, _, at_center, _ = face_enclosures(centers, det, derivatives)
        d, de, _, gradient = face_enclosures(box, det, derivatives)
        action = at_center
        for g, (a,b) in zip(gradient, box):
            action += g * arb(0, (b-a)/2)
        previous = box.copy()
        for i, (a,b) in enumerate(previous):
            derivative = gradient[i] - action * de[i]/(4*d)
            if derivative > 0:
                box[i] = (b,b)
            elif derivative < 0:
                box[i] = (a,a)
        assert box != previous, "a maximizing coordinate needs subdivision"
    d, _, action, _ = face_enclosures(box, det, derivatives)
    return 384*action/d.sqrt().sqrt()


def main():
    ctx.prec = 192
    det = determinant_polynomial()
    derivatives = [det.derivative(i) for i in range(10)]
    dmin, dmax = determinant_extrema(det, derivatives)
    m, M = [arb(d).sqrt().sqrt() for d in (dmin, dmax)]
    scale = fmpq(823935,100000)
    assert 16*M < scale
    local = local_maximum(det, derivatives)
    assert local < fmpq(-13692838508, 10**8)
    global_upper = local*2*(m*M).sqrt()/(m+M)
    b = fmpq(136914,1000)
    assert global_upper < -b

    norm = (120*6**5*arb(fmpq(7,8))**9)**256*(15*arb.pi())**640/factorial(1920)
    volume = arb(1)
    for (lo,hi), count in zip(GAPS,(1024,1536,1024,256)):
        volume *= (arb(hi-lo)/(2000*scale))**count
    p = volume/norm
    assert p > 2*arb(10)**-7013
    a = 5*arb(6144).sqrt()/4
    chebyshev = (1/p-1).acosh()/(b/a).acosh()
    monomial = (1/p).log()/(b/a).log()
    assert 18682 < chebyshev < 18683
    assert 48259 < monomial < 48260
    print("global shape upper bound", global_upper)
    print("minus log10 probability bound", -p.log()/arb(10).log())
    print("Chebyshev and monomial degree thresholds", chebyshev, monomial)


if __name__ == "__main__":
    main()
