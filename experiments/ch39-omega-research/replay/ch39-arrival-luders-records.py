"""Sharp number of channel-diagonal records for the paired Fourier arrival clock."""
from collections import Counter
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
from runpy import run_path

from sympy import Poly, cyclotomic_poly, symbols

root = Path(__file__).parent
clock = run_path(str(root/'ch39-charged-arrival-clock.py'))
shell = clock['shell']


def local_support():
    squared = shell['squared']
    adjacency = {p:{q:w for q,w in row if squared(p,q)>=108}
                 for p,row in shell['graph'](3,56).items()}
    pole = (51,51,51)
    adjacency[pole] = {p:0 for p in adjacency if 108<=squared(p,pole)<=504}
    for p in adjacency[pole]:
        adjacency[p][pole]=0
    reach = {}
    for i, s in enumerate(clock['sites']):
        layer = {(clock['orbit'](s)[0],0)}
        for _ in range(5):
            layer = {(v,k+w) for p,k in layer for v,w in adjacency[p].items()}
        for j, t in enumerate(clock['sites']):
            reach[i,j] = {d for d,v in enumerate(clock['orbit'](t))
                          if all((v,k) in layer for k in (-1,0,1))}
    coords = Counter()
    for o,h,a in product(range(3),range(16),range(13)):
        i = clock['position'](h,a)
        # Second input pair uses direction +1 for the partner's long digits.
        shift = int(12<=h<16 and a<4)
        coords[i,(o+shift)%3] += 1
        for b in range(13):
            j = clock['position'](h^8,b)
            target_shift = int(12<=(h^8)<16 and b<4)
            assert (target_shift-shift)%3 in reach[i,j]
    assert all(coords[i,o]==n for i,n in enumerate(clock['per_site']) for o in range(3))


def verify():
    local_support()
    t = symbols('t')
    modulus = Poly(cyclotomic_poly(13,t),t)
    for d in range(1,13):
        assert Poly(sum(t**((a*d)%13) for a in range(13)),t).rem(modulus).is_zero
        assert not Poly(sum(t**((a*d)%13) for a in range(4)),t).rem(modulus).is_zero
    p, q = Q(4,13), Q(9,13)
    assert p*q*q-p*p*q == Q(180,2197)
    print('48 packets of 13, all 13-by-13 packet transitions have five-edge support.')
    print('Every off-diagonal digit pair has nonzero long-projector Gram entry.')
    print('13 digit records suffice. Purely binary retracing excess at least', p*q*(q-p))


if __name__=='__main__':
    verify()
