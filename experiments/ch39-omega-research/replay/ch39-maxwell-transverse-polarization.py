"""Finite-difference Maxwell symbol for one transverse polarization in the two-pair family."""
from itertools import product
import sympy as sp


def transverse():
    a,b = sp.symbols('a b', nonzero=True)
    gram = sp.zeros(3)
    def mode(s):
        return sp.Matrix([int(len(s & {2,3}) == k)*(int(0 in s)-int(1 in s)) for k in range(3)])
    for colors in product(range(3), repeat=4):
        s, t = ({i for i,c in enumerate(colors) if c == j} for j in (1,2))
        if not s or not t:
            continue
        delta = mode(s)+a**len(s & {0,1})*b**len(s & {2,3})*mode(t)-mode(s|t)
        gram += delta.xreplace({a:1/a,b:1/b})*delta.T
    return a,b,sp.factor(gram.det())


def verify():
    a,b,det = transverse()
    u,v,T,V = sp.symbols('u v T V', nonzero=True)
    conjugate = det.xreplace({a:b,b:a})
    change = {a:u*v,b:u/v}
    even = sp.expand(((det+conjugate)/32).subs(change,simultaneous=True))
    odd = ((det-conjugate)/32).subs(change,simultaneous=True)
    def pair(n,x):
        return 1 if n==0 else 2*sp.chebyshevt(n,x/2)
    symmetric = 0
    for term in sp.Add.make_args(even):
        powers = term.as_powers_dict()
        i,j = powers[u],powers[v]
        if i>=0 and j>=0:
            symmetric += term/(u**i*v**j)*pair(i,T)*pair(j,V)
    assert sp.cancel(symmetric.subs({T:u+1/u,V:v+1/v})-even)==0
    antisymmetric = (2*T**4-20*T**2+88+(T**3-14*T)*V)
    assert sp.cancel(odd+(u-1/u)*(v-1/v)*antisymmetric.subs({T:u+1/u,V:v+1/v})/2)==0
    on_cross = sp.factor(symmetric.subs(V,12/T-2*T))
    assert sp.factor(on_cross-(T*T-4)**2*(T**4-11*T*T+108)/T**2)==0
    assert sp.factor(antisymmetric.subs(V,12/T-2*T)-20*(T*T-4))==0
    # A positive real part excludes either transverse Maxwell zero on this branch.
    assert 11**2 < 4*108
    print('Even Maxwell determinant on the cross graviton branch:', on_cross)


if __name__ == '__main__':
    verify()
