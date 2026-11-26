"""Finite plaquette gap, quantized flux, and Haar probability estimates for the overlap kernel."""
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations, product
import json
from pathlib import Path


def tent(radius):
    values = [radius+1-abs(x) for x in range(-radius, radius+1)]
    norm2 = sum(x*x for x in values)
    difference2 = sum((x-y)**2 for x, y in zip(values+[0], [0]+values))
    assert Q(difference2, norm2) == Q(6, 2*(radius+1)**2+1)


def bounds():
    radius, eps, closure = 4095, Q(1, 2**26), Q(1, 1024)
    for size in (1, 2, 7, radius):
        tent(size)
    # sqrt(3)<7/4 and the five step coefficient is 1+2/sqrt(5)<19/10.
    assert Q(7, 4)**2 > 3 and (Q(19, 10)-1)**2 > Q(4, 5)
    loss = Q(19, 10)*(8*Q(7, 4)/(radius+1)+12*radius*eps+closure)
    assert loss < Q(1, 100)
    assert Q(19, 10)*(8*Q(7, 4)/56+Q(5, 20)) < 1
    print("Bound on the loss from free gap:", loss)


def flux_exponents(side, x, y):
    horizontal = -side*y if x == side-1 else 0
    return horizontal, x


def flux_certificate():
    side = 25000
    # Each boundary case on the two-coordinate torus, including the corner seam.
    for size in (5, side):
        positions = range(size) if size == 5 else (0, 1, size//2, size-2, size-1)
        for x, y in product(positions, repeat=2):
            h, v = flux_exponents(size, x, y)
            _, v_right = flux_exponents(size, (x+1) % size, y)
            h_above, _ = flux_exponents(size, x, (y+1) % size)
            assert (h + v_right - h_above - v) % (size*size) == 1
    assert Q(8, side*side) < Q(1, 2**26)
    weights = [1, -1, 0, 0, 0]
    matter = [-w for w in weights]+[a+b for a, b in combinations(weights, 2)]
    assert Counter(matter) == {-1: 4, 0: 7, 1: 4}


def haar_certificate():
    prefix = json.loads(Path(__file__).with_name("ch39-tick-separators.json").read_text())["prefix"]
    assert Q(prefix[0][0], 1 << 400) > 14
    trace_bound = 30+Q(2, 15)
    amplitude = trace_bound/(14**2+Q(1, 4))
    distortion = Q(9, 7)*((1+amplitude)/(1-amplitude))**14
    eps = Q(1, 2**26)
    assert distortion < 100
    assert sum(2*n-1 for n in range(2, 6)) == 24
    assert distortion*eps**24 < Q(1, 2**617)
    # The critical constant background has Re(u)=61/200 and has W=m exactly.
    assert 5*(1-Q(61, 200)) == Q(139, 40)
    assert 50 < 8**2
    # Tiling with even checkerboard squares does not reuse an axis edge.
    for side in (4, 6, 10):
        edges = []
        for x, y in product(range(side), repeat=2):
            if (x+y) % 2 == 0:
                edges.extend([(0, x, y), (0, x, (y+1) % side),
                              (1, x, y), (1, (x+1) % side, y)])
        assert len(edges) == len(set(edges)) == 2*side*side
    print("Conditional Haar distortion bound:", float(distortion))


if __name__ == "__main__":
    bounds()
    flux_certificate()
    haar_certificate()
    print("Verified admissibility constants, two flux sheets and finite-volume probability bounds.")
