"""Exact cross polarization on two equal pairs of lattice momentum components."""
from collections import defaultdict
from itertools import permutations
from pathlib import Path
from runpy import run_path
import sympy as sp

base = run_path(str(Path(__file__).with_name('ch39-regge-flat-cell.py')))
EDGES = base['EDGES']


def cross_symbol():
    local = sum((m for _,_,m in base['canonical']()),sp.zeros(10))
    mode = {5:1, 6:-1, 9:-1, 10:1}
    stencil = defaultdict(int)
    # All four edges supporting the cross mode have the same midpoint Fourier offset.
    for order in permutations(range(4)):
        vertices = [sum(1 << i for i in order[:j]) for j in range(5)]
        for (i,j), value in local.todok().items():
            a,b = [vertices[k] for k in EDGES[i]]
            c,d = [vertices[k] for k in EDGES[j]]
            if b-a not in mode or d-c not in mode:
                continue
            t = c.bit_count()-a.bit_count()
            s = (sum(((c>>k)&1)-((a>>k)&1) for k in (0,1))
                 -sum(((c>>k)&1)-((a>>k)&1) for k in (2,3)))
            stencil[t,s] += mode[b-a]*mode[d-c]*value/4
    T,V = sp.symbols('T V')
    def cosines(n, x):
        return 1 if n == 0 else 2*sp.chebyshevt(n,x/2)
    assert all(val==stencil[-t,s]==stencil[t,-s] for (t,s),val in list(stencil.items()))
    return sp.factor(sum(val*cosines(t,T)*cosines(s,V) for (t,s),val in stencil.items() if t>=0 and s>=0))


def verify():
    T,V = sp.symbols('T V')
    spin2 = cross_symbol()
    assert sp.expand(spin2-(2*T**2+T*V-12)/20) == 0
    scalar = 32-(T**2-2)*(T+V)**2
    residual = sp.factor(scalar.subs(V,12/T-2*T))
    assert sp.factor(residual-(T*T-4)**2*(18-T*T)/T**2) == 0
    # Formal small spatial momentum: E^2=s^2/5+c*s^4.
    s,c = sp.symbols('s c')
    time = 2+s**2/5+(c+sp.Rational(1,300))*s**4
    space = 2-s**2+s**4/12
    for expression, expected in [(spin2,-sp.Rational(1,125)),(scalar,sp.Rational(1,50))]:
        series = expression.subs({T:time,V:space}).expand()
        assert series.coeff(s,0)==series.coeff(s,2)==0
        assert sp.solve(series.coeff(s,4),c)==[expected]
    print('Cross graviton symbol:', spin2, 'scalar evaluated on its positive transfer branch:', residual)


if __name__ == '__main__':
    verify()
