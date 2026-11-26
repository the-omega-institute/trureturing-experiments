"""Exact vertex conformal Hessian, a negative minor bound, and the period-four fibers."""
from collections import Counter
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import fmpz_mat

symbol = run_path(str(Path(__file__).with_name("ch39-regge-fierz-pauli.py")))


def step(mask):
    return tuple((mask >> i) & 1 for i in range(4))


def vertex_coefficients(midpoint):
    pulled, anchored, bound = Counter(), Counter(), Counter()
    for (i,j,delta), value in midpoint.items():
        a,b = step(i+1),step(j+1)
        # Displace the two edge midpoints back to their respective starts.
        twice = tuple(d+x-y for d,x,y in zip(delta,a,b))
        assert all(v % 2 == 0 for v in twice)
        start = tuple(v//2 for v in twice)
        anchored[i,j,start] += value
        bound[i+1,j+1] += abs(value)
        weight = value*sum(a)*(5-sum(a))*sum(b)*(5-sum(b))
        for left,right in product((0,1), repeat=2):
            offset = tuple(d-left*x+right*y for d,x,y in zip(start,a,b))
            pulled[offset] += weight
    expected = Counter({(0,)*4: 63360})
    for mask in range(1,16):
        v = step(mask)
        for sign in (-1,1):
            expected[tuple(sign*x for x in v)] = -4608 if sum(v) in (1,4) else -864
    assert pulled == expected  # H was scaled by 240, the two q factors by 25.
    minor = (3,5,9,14)
    for a in minor:
        diagonal = {(d,v) for (i,j,d),v in anchored.items() if i == j == a-1 and v}
        assert diagonal == {((0,)*4,-96)}
        assert sum(bound[a,b] for b in minor if b != a) == 72
    return anchored


def sign_changes(coefficients):
    signs = [c > 0 for c in coefficients if c]
    return sum(a != b for a,b in zip(signs,signs[1:]))


def inertia(anchored, momentum):
    # Endpoint phases on a period-four box are Gaussian integers.
    block = [[0]*30 for _ in range(30)]
    for (i,j,offset),v in anchored.items():
        phase = sum(x*k for x,k in zip(offset,momentum)) % 4
        re,im = ((1,0),(0,1),(-1,0),(0,-1))[phase]
        block[i][j] += re*v
        block[i+15][j+15] += re*v
        block[i+15][j] += im*v
        block[i][j+15] -= im*v
    real = fmpz_mat(block)
    assert real == real.transpose()
    poly = list(real.charpoly())  # ascending powers
    zero = next(i for i,c in enumerate(poly) if c)
    pos = sign_changes(poly)
    neg = sign_changes([(-1)**i*c for i,c in enumerate(poly)])
    assert pos+neg+zero == 30 and all(x % 2 == 0 for x in (pos,neg,zero))
    return (pos//2, neg//2, zero//2)


def verify():
    anchored = vertex_coefficients(symbol["coefficients"]())
    signatures = Counter(inertia(anchored,k) for k in product(range(4),repeat=4))
    print("Period-four momenta by (positive,negative,zero) inertia:", signatures)
    assert signatures == {(0,5,10):1, (1,10,4):255}


if __name__ == "__main__":
    verify()
