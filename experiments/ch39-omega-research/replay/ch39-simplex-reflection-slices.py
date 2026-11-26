"""Primitive reflection slices of the four-simplex lattice and finite two-layer positivity."""
from itertools import permutations, product
import sympy as sp


def slice_checks():
    step = [sp.eye(4)[:,j] for j in range(4)]+[-sp.ones(4,1)]
    g = sp.eye(4)-sp.ones(4)/5
    assert sum(step[1:],step[0]) == sp.zeros(4,1)
    normals = set()
    for i in range(5):
        for j in range(i+1,5):
            d = step[i]-step[j]
            h = d.T*g
            assert h*d == sp.Matrix([[2]]) and (d.T*g*d)[0]==2
            reflection = sp.eye(4)-d*h
            assert reflection.T*g*reflection==g and reflection.det()==-1
            assert all(reflection*v==step[j if k==i else i if k==j else k] for k,v in enumerate(step))
            values = [(h*v)[0] for v in step]
            for order in permutations(values):
                path = [sum(order[:k]) for k in range(5)]
                assert max(path)-min(path)==1
            normals.add(tuple(h))
    assert len(normals)==10
    for u in product(range(-2,3), repeat=4):
        values = (*u,-sum(u))
        if sum(max(0,x) for x in values)==1:
            assert u in normals or tuple(-x for x in u) in normals


def scalar_symbol():
    s,t,u,v = sp.symbols("s t u v", nonzero=True)
    phases = (s*t,s/t,u,v)
    def conjugate(expr):
        return expr.xreplace({x:1/x for x in (s,t,u,v)})
    half = s*(1+u)*(1+v)
    expected = 32-(half+conjugate(half))*(s+1/s+t+1/t)
    action = 0
    for bits in product(range(2),repeat=4):
        monomial = sp.prod(p**b for p,b in zip(phases,bits))
        action += (1-monomial)*(1-1/monomial)
    assert sp.expand(action-expected)==0


def transfer_example():
    # Each boundary has a five-state variable. This positive pointwise face
    # weight fails the separate one-face Fourier test, while gluing gives a Gram matrix.
    values = [1,1,0,0,1]
    layer = sp.Matrix(5,5,lambda i,j: values[(i-j)%5])
    assert any(x<0 for x in layer.eigenvals())
    glued = layer*layer.T
    assert all(x>0 for x in glued.eigenvals())
    assert glued[0,0]==3 and glued[0,1]==2


if __name__ == "__main__":
    slice_checks()
    scalar_symbol()
    transfer_example()
    print("Ten primitive slices, the scalar Laurent identity, and finite one- versus two-layer positivity verified.")
