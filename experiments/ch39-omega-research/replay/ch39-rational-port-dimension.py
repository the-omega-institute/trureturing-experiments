"""Minimum rational channel dimension and a one-register realization of gain 624."""
from collections import Counter
from itertools import product
from math import isqrt

from sympy import I, Matrix, cancel, eye, symbols


def transfer():
    z = symbols("z")
    v = Matrix([1,24,4,4,4])
    H = 2*v*v.T/(v.dot(v))-eye(5)
    assert H.T*H == eye(5)
    a, b, c, D = H[0,0], H[:1,1:], H[1:,:1], H[1:,1:]
    result = a+(z*b*(eye(4)-z*D).inv()*c)[0]
    target = ((625*z)-623)/(625-623*z)
    assert cancel(result-target) == 0
    assert cancel((result-1)/(result+1)-624*(z-1)/(z+1)) == 0
    complex_column = Matrix([1,24+4*I,4+4*I])
    norm = (complex_column.conjugate().T*complex_column)[0].expand()
    assert norm == 625
    unitary = 2*complex_column*complex_column.conjugate().T/norm-eye(3)
    assert (unitary*unitary.conjugate().T-eye(3)).applyfunc(cancel) == Matrix.zeros(3)
    assert 9*8**2+3*4**2 == 624
    assert 39%8 not in {sum(t*t for t in triple)%8 for triple in product(range(8),repeat=3)}


def register():
    pole = (17,17,17)
    shell = {p for p in product(range(30),repeat=3) if sum(x*x for x in p)==867}
    assert len(shell) == 19
    near = (13,13,23)
    def rotate(p):
        return p[1],p[2],p[0]
    cycle = (near,rotate(near),rotate(rotate(near)))
    states = [(p,j) for p in sorted(shell) for j in range(3 if p in cycle else 1)]
    inner = [s for s in states if s != (pole,0)]
    def lift(s):
        p,j = s
        return rotate(p), (j+int(p==cycle[2]))%3
    weights = {s:6 if len(set(s[0]))<3 else 4 for s in inner}
    assert all(weights[lift(s)] == w for s,w in weights.items())
    assert Counter(weights.values()) == {6:12,4:12}
    assert sum(w*w for w in weights.values()) == 624
    def mirror(s):
        p,j = s
        if p in cycle:
            phase = (-(cycle.index(p)+3*j))%9
            return cycle[phase%3],phase//3
        return (p[1],p[0],p[2]),j
    assert all(mirror(mirror(s)) == s and lift(mirror(lift(s))) == mirror(s) for s in inner)
    assert all(weights[mirror(s)] == weights[s] for s in inner)
    solutions = [(i,j) for i in range(1,isqrt(52)+1) for j in range(i,isqrt(52)+1)
                 if i*i+j*j==52]
    assert solutions == [(4,6)]
    v = Matrix([1,*weights.values()])
    H = 2*v*v.T/625-eye(25)
    assert H*H == eye(25)
    # The perturbed phase is shared by a full class; the permutation fixes its port weight.
    marked_one = sum(w*w for s,w in weights.items() if len(set(s[0]))<3)
    marked_two = sum(1 for a,b in product(states,states) if a in weights and len(set(a[0]))<3)
    assert (marked_one,marked_two)==(432,300)
    return len(inner), marked_one, marked_two


if __name__ == "__main__":
    transfer()
    print("Minimum real rational internal dimension 4, shell internal states and marked gains",register())
