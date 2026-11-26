"""Stabilizers and CP capacity of the shell with an optional tetrahedron orientation."""
from collections import Counter
from fractions import Fraction as F
from itertools import permutations, product

import sympy as sp


def sign(p):
    return (-1)**sum(u>v for i,u in enumerate(p) for v in p[i+1:])


def fields(p,chi):
    h = tuple(p[i]+F(chi,128)*(p[(i+1)%3]-p[(i+2)%3]) for i in range(3))
    delta = (h[0]-h[1])*(h[1]-h[2])*(h[2]-h[0])
    t = delta/2**16
    assert min(h)>0 and abs(t)<1
    return h,t


def verify():
    shell = [p for p in product(range(30),repeat=3) if sum(v*v for v in p)==867]
    perms = list(permutations(range(3)))
    count = Counter()
    mean_det = F(0)
    for p,chi in product(shell,(0,-1,1)):
        h,t = fields(p,chi)
        stab = 0
        for perm in perms:
            image = tuple(p[i] for i in perm)
            orientation = sign(perm)*chi
            hh,tt = fields(image,orientation)
            assert hh == tuple(h[i] for i in perm) and tt==sign(perm)*t
            stab += (image,orientation)==(p,chi)
        assert (len(set(h))==3) == (stab==1)
        count[chi,stab] += 1
        mean_det += 2*(t*2**16)*(3*t-t**3)/57
    assert count == {(0,1):12,(0,2):6,(0,6):1,(-1,1):18,(-1,3):1,(1,1):18,(1,3):1}
    assert mean_det>0
    x,y,z,t = sp.symbols("x y z t", real=True)
    U = sp.diag(x,y,z)
    A = sp.Matrix([[0,1,-1],[-1,0,1],[1,-1,0]])
    D = 4*sp.eye(3)+sp.ones(3)+sp.I*t*A
    assert sp.simplify((U*D-D*U).det()-2*sp.I*(x-y)*(y-z)*(z-x)*(3*t-t**3)) == 0
    print("Allowed nonzero CP determinant: 12+18+18 of 57. Signed determinant mean:",mean_det)


if __name__ == "__main__":
    verify()
