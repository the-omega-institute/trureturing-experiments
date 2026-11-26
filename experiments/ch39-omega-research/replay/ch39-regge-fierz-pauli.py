"""Exact midpoint Bloch symbol through second order and the linearized Einstein kinetic tensor."""
from collections import defaultdict
from itertools import combinations, permutations
from pathlib import Path
from runpy import run_path
import sympy as sp

base = run_path(str(Path(__file__).with_name('ch39-regge-flat-cell.py')))
EDGES = base['EDGES']
METRIC = [(i,j) for i in range(4) for j in range(i,4)]


def coefficients():
    local = sum((v for _, _, v in base['canonical']()), sp.zeros(10))
    assert local == local.T
    coeff = defaultdict(int)
    for perm in permutations(range(4)):
        vertices = [sum(1 << i for i in perm[:k]) for k in range(5)]
        for (i,j), val in local.todok().items():
            v,w = [vertices[p] for p in EDGES[i]]
            x,y = [vertices[p] for p in EDGES[j]]
            delta = tuple(((x>>k)&1)+((y>>k)&1)-((v>>k)&1)-((w>>k)&1) for k in range(4))
            a,b = w-v-1,y-x-1
            assert (240*val).is_Integer
            coeff[a,b,delta] += int(240*val)
    coeff = {k:v for k,v in coeff.items() if v}
    assert all(v == coeff[a,b,tuple(-n for n in delta)] for (a,b,delta),v in coeff.items())
    # In midpoint convention the gauge vector for vertex motion is 2i sin(k.S/2) S.
    gauge = defaultdict(int)
    for (a,b,delta), value in coeff.items():
        step = tuple((b+1)>>k & 1 for k in range(4))
        for sign in (-1,1):
            exponent = tuple(d+sign*s for d,s in zip(delta,step))
            for k, bit in enumerate(step):
                gauge[a,k,exponent] += sign*value*bit
    assert not any(gauge.values())
    return coeff


def continuum(k):
    inverse = sp.eye(4)+sp.ones(4)
    k = sp.Matrix(k)
    raised = inverse*k
    kk = (k.T*raised)[0]
    basis = []
    for i,j in METRIC:
        h = sp.zeros(4)
        h[i,j] = h[j,i] = 1
        basis.append(h)
    def fp(h, v):
        th, tv = (inverse*h).trace(), (inverse*v).trace()
        return (kk*((inverse*h*inverse*v).trace()-th*tv)-2*(raised.T*h*inverse*v*raised)[0]
                +th*(raised.T*v*raised)[0]+tv*(raised.T*h*raised)[0])
    return sp.Matrix([[fp(h,v) for v in basis] for h in basis])


def verify():
    coeff = coefficients()
    modes = sp.Matrix(15,10,lambda a,b:(1+int(METRIC[b][0]!=METRIC[b][1]))*
                      int(all((a+1)&(1<<j) for j in METRIC[b])))
    ratios = set()
    for i,j in METRIC:
        k = [int(d==i or d==j) for d in range(4)]
        quadratic = sp.zeros(15)
        for (a,b,delta), v in coeff.items():
            quadratic[a,b] -= v*sum(s*x for s,x in zip(k,delta))**2
        pulled = modes.T*quadratic*modes/(8*240)
        fp = continuum(k)
        pivot = next(pos for pos in fp.todok())
        ratio = pulled[pivot]/fp[pivot]
        assert pulled == ratio*fp
        ratios.add(ratio)
    assert ratios == {-sp.Rational(1,20)}
    print('Even midpoint Laurent entries, exact gauge identities, projected Fierz-Pauli factor:', len(coeff), ratios)


if __name__ == '__main__':
    verify()
