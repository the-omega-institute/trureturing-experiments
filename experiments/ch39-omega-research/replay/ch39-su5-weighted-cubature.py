"""Exact fourth Haar moment for ten SU(2) factors and bounded positive face cubature."""
from collections import defaultdict
from fractions import Fraction as F
from functools import cache
from itertools import product
from math import factorial
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

full = run_path(str(Path(__file__).with_name("ch39-full-color-obstruction.py")))


def entries(axis, left, right, conjugate):
    # Each option carries its output index, phase powers for a and b, and sign.
    if axis == left:
        options = [(left,1,0,1), (right,0,1,1)]
    elif axis == right:
        options = [(left,0,-1,-1), (right,-1,0,1)]
    else:
        options = [(axis,0,0,1)]
    if conjugate:
        return [(j,-a,-b,sign) for j,a,b,sign in options]
    return options


@cache
def column(indices, left, right, m):
    options = [entries(axis,left,right,k>=2) for k,axis in enumerate(indices)]
    result = defaultdict(F)
    for terms in product(*options):
        out, a, b, signs = zip(*terms)
        if sum(a) or sum(b):
            continue
        i, j = sum(map(abs,a))//2, sum(map(abs,b))//2
        sign = (-1)**signs.count(-1)
        result[out] += sign*F(m*factorial(m+i-1)*factorial(j),factorial(m+i+j))
    return dict(result)


def fourth_moment():
    # Act on a tensor vector from right to left, peeling the smallest sphere first.
    factors = [(left,right,5-right) for left in reversed(range(4)) for right in range(left+1,5)]
    pairs = list(product(range(5), repeat=2))
    for source in product(range(5), repeat=4):
        state = {source:F(1)}
        for left,right,m in factors:
            target = defaultdict(F)
            for key,amplitude in state.items():
                for out,value in column(key,left,right,m).items():
                    target[out] += value*amplitude
            state = {key:value for key,value in target.items() if value}
        # Project onto the two invariant contractions with their exact inverse Gram matrix.
        i,j,k,l = source
        c1,c2 = int((i,j)==(k,l)), int((i,j)==(l,k))
        expected = defaultdict(F)
        for x,y in pairs:
            expected[x,y,x,y] += F(5*c1-c2,120)
            expected[x,y,y,x] += F(5*c2-c1,120)
        assert state == {key:value for key,value in expected.items() if value}


def face_bounds():
    ctx.prec = 160
    q = full["q"]
    strength = q(1)/(q(11134)**2 + q(F(1,4)))
    trace_bound = q(F(452,15))
    max_u = 2*strength*trace_bound/(1+strength*trace_bound)
    max_power = q(30484477)
    assert max_power*max_u < 16 and (1-max_u)**max_power > q(F(1,2**23))
    error = F((2**23)*16**97, factorial(97))
    assert 8*10**6*error < F(1,10**20)
    assert 4+2*14*96 == 2692
    assert (2692//2+2)//2 == 674


if __name__ == "__main__":
    fourth_moment()
    face_bounds()
    print("All 625 quartic tensor columns and the amplified finite-volume cubature margins passed.")
