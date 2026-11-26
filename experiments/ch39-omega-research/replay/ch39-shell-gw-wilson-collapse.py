"""Exact singular-gain asymptotics near the shortest-arc action minima."""
from itertools import product
import sympy as sp


def leading_columns(n):
    # A cyclic choice of the two moving axes at each fixed axis is equivalent
    # to the other orientation since both words and their adjoints are present.
    s=sp.Symbol("s",real=True)
    eye=sp.eye(3)
    directions=[]
    for k in range(n):
        i,j=(k+1)%3,(k+2)%3
        d=sp.diag(*[1 if a==k else -1 for a in range(3)])
        a=sp.zeros(3)
        a[i,j],a[j,i]=-1,1
        directions += [d+s*a,d-s*a]
    words=[eye]+directions+[a*b for a,b in product(directions,repeat=2)]
    columns=[]
    degrees=[]
    for i,j in product(range(3),repeat=2):
        polynomials=[sp.Poly(word[i,j],s) for word in words]
        degree=next(d for d in range(3) if any(p.nth(d)!=0 for p in polynomials))
        columns.append(sp.Matrix([p.nth(degree) for p in polynomials]))
        degrees.append(degree)
    return columns,degrees


def residual_gram(n):
    cols,deg=leading_columns(n)
    order=max(deg)
    low=sp.Matrix.hstack(*(v for v,d in zip(cols,deg) if d<order))
    high=sp.Matrix.hstack(*(v for v,d in zip(cols,deg) if d==order))
    gram=low.T*low
    assert gram.det()>0
    return order,high.T*high-high.T*low*gram.inv()*low.T*high


if __name__=="__main__":
    for n in (2,3):
        degree,gram=residual_gram(n)
        assert (degree,gram.eigenvals()) == ((2,{4:2}) if n==2 else (1,{32:3,20:3}))
        print(n,degree,gram.eigenvals())
