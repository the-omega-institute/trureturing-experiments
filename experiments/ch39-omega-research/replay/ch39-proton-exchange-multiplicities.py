"""Minimum extra triplets needed to kill successive tree-level momentum coefficients."""
from fractions import Fraction as Q
from itertools import combinations, product
from pathlib import Path
from runpy import run_path

root = Path(__file__).resolve().parent
old = run_path(str(root/'ch39-eight-time-pole-control.py'))
B, mul, rounded = (old[key] for key in ('B','mul','rounded'))
ids = [1,2,3,4,8,9,51,52]
nodes,weights = old['quadrature'](ids,Q(0))
x = dict(zip(ids,[rounded(v/nodes[3]) for v in nodes],strict=True))
legs = dict(zip(ids,[rounded(w*v) for w,v in zip(weights,nodes)],strict=True))


def absolute(v):
    if v.lo >= 0: return v
    if v.hi <= 0: return B(-v.hi,-v.lo)
    return B(0,max(-v.lo,v.hi))


def product_at(i, roots):
    if i in roots: return B(0)
    return mul(*(x[i]-x[j] for j in roots))


def check_necessary(k, exclusions):
    for i,roots in exclusions.items():
        assert i not in roots and len(set(roots))==len(roots)<k
        sizes = {j:absolute(legs[j]*product_at(j,roots)) for j in ids}
        assert sizes[i].lo>sum(v.hi for j,v in sizes.items() if i!=j)


def poly_with_roots(roots):
    coeff = [B(1)]
    for i in roots:
        previous = coeff
        coeff = [B(0)]*(len(previous)+1)
        for j,v in enumerate(previous):
            coeff[j] = rounded(coeff[j]-x[i]*v)
            coeff[j+1] = rounded(coeff[j+1]+v)
    return coeff


def sign(v):
    if v.lo>0: return 1
    assert v.hi<0
    return -1


def strictly_positive_polynomial(coeff):
    # Sturm divisions use the exact algebraic cancellation of each leading term.
    assert len(coeff)%2==1 and sign(coeff[-1])==1
    a = coeff
    b = [rounded(B(i)*v) for i,v in enumerate(a) if i]
    seq = [a,b]
    while len(b)>1:
        r = list(a)
        while len(r)>=len(b):
            shift = len(r)-len(b)
            factor = rounded(r[-1]/b[-1])
            for j in range(len(b)-1):
                r[shift+j] = rounded(r[shift+j]-factor*b[j])
            r.pop()
        a,b = b,[B(-v.hi,-v.lo) for v in r]
        sign(b[-1])
        seq.append(b)
    right = [sign(p[-1]) for p in seq]
    left = [s*(-1)**(len(p)-1) for s,p in zip(right,seq)]
    variations = lambda signs: sum(a!=b for a,b in zip(signs,signs[1:]))
    assert variations(left)==variations(right)


def construct_complex(k, simple, extra):
    # A nonnegative polynomial U=|P|^2 fixes the moduli of all simple residues.
    d = {i:legs[i]*absolute(product_at(i,[j for j in ids if j!=i])) for i in ids}
    squared = {i:rounded((v/d[52])**2) for i,v in d.items()}
    coeff = [B(0)]*(len(simple)+bool(extra))
    for i in simple:
        others = [j for j in simple if j!=i]
        scale = rounded(squared[i]/product_at(i,others))
        for j,v in enumerate(poly_with_roots(others)):
            coeff[j] = rounded(coeff[j]+scale*v)
    if extra:
        for j,v in enumerate(poly_with_roots(simple)):
            coeff[j] = rounded(coeff[j]+B(extra)*v)
    assert len(coeff)-1<=2*(len(ids)-k-1)
    strictly_positive_polynomial(coeff)
    for i in set(ids)-set(simple):
        val = B(0)
        for c in reversed(coeff): val = rounded(val*x[i]+c)
        assert val.hi<squared[i].lo


def classify_complex():
    forced = {
        1:{1:[]}, 2:{1:[]},
        3:{1:[],2:[1,8]}, 4:{1:[],2:[1,8]},
        5:{1:[],2:[1,8],3:[1,2,8,9]},
        6:{1:[],2:[1,8],3:[1,2,8,9],4:[1,2,3,8,52],8:[1,2,3,4,52],52:[1,2,3,4,8]},
        7:{i:[j for j in ids if j not in (i,51)] for i in ids if i!=51}}
    for k,mapping in forced.items(): check_necessary(k,mapping)
    b = [legs[j]*absolute(x[1]-x[j]) for j in ids[1:]]
    for i in range(len(b)): assert b[i].hi<sum(v.lo for j,v in enumerate(b) if i!=j)
    construct_complex(4,[3,4,8,9,51,52],33000000)
    construct_complex(5,[4,8,9,51,52],0)
    # At orders six and seven the forced count already equals k, and the
    # vertex of the real cube section gives k additional copies as an upper bound.
    print('Minimum additional copies for k=1,...,8, complex couplings: 1 1 2 2 3 6 7 8')


def classify_real():
    for k in range(1,9):
        for doubled in combinations(ids,k-1):
            simple = [i for i in ids if i not in doubled]
            values = [legs[i]*product_at(i,doubled) for i in simple]
            for signs in product((-1,1),repeat=len(values)-1):
                test = values[0]+sum((B(s)*v for s,v in zip(signs,values[1:])),B(0))
                assert test.lo>0 or test.hi<0
    print('Minimum additional copies with real couplings: 1 2 3 4 5 6 7 8')


if __name__=='__main__':
    classify_complex()
    classify_real()
