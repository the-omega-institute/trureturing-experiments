"""Axis potential of the fully split scalar and Dirac families on the original arcs."""
from fractions import Fraction
from pathlib import Path
import runpy
from types import FunctionType
from flint import arb

root = Path(__file__).resolve().parent
mass = runpy.run_path(str(root/'ch39-su5-wilson-mass.py'))
u8 = runpy.run_path(str(root/'ch39-unit-pentagon.py'))
m,q = u8['model'],u8['q']
phase = mass['phase']
net = {2:8,3:6,7:1,1:3}
sectors = {}
for p in net:
    shell = dict(m['shell'])
    shell['edges'] = [(i,j,d,z**p) for i,j,d,z in shell['edges']]
    f = m['schur_terms']
    sectors[p] = FunctionType(f.__code__,dict(f.__globals__,shell=shell),f.__name__,f.__defaults__,f.__closure__)

def slopes(kappa,cosines):
    degree=20
    cutoff = 8192
    first=Fraction(1,1024)
    weights = [(q(0).union(q(first)),q(first**2/2))]
    for power in range(-10,13):
        half=q(Fraction(2)**power)/2
        for j in range(degree):
            abscissa, w = arb.legendre_p_root(degree,j,weight=True)
            k=half*(3+abscissa)
            weights.append((k,k*half*w))
    totals = [q(0) for c in cosines]
    error = 4*q(6000)*cutoff/2**(2*degree) + m['tail_bound'](cutoff)
    for p,dim in net.items():
        delta = kappa*(1-(phase**p).real)
        for k,w in weights:
            lo,hi = u8['node_enclosure'](k+delta)
            terms = sectors[p](lo,hi)
            for j,c in enumerate(cosines):
                totals[j] += dim*w*m['log_derivative'](terms,c)
    return [t+arb(0,(sum(net.values())*error).upper()) for t in totals]


def static_gap():
    flow = mass['flow']
    for ch in (1,2,6):
        matrix = mass['m']['readers'][ch].__globals__['low_matrix']
        for j in range(64):
            half = arb.pi()/128
            h = matrix(mass['m']['u'],mass['m']['v'],(2*j+1)*half)
            assert flow['inertia'](h,q('.25'))==0
    assert q('.25')-4*(arb.pi()/128+q('.00001'))-flow['high_error']>q('.08')
    # The original p=1 static floor follows from the earlier inner triangle certificate.

if __name__=='__main__':
    static_gap()
    cs = [q(Fraction(i-8,8)).union(q(Fraction(i-7,8))) for i in range(16)]
    kappa = q('2.738').union(q('2.755'))
    all_slopes = slopes(kappa,cs)
    lower = min(v.lower() for v in all_slopes)
    mesh = runpy.run_path(str(root/'ch39-gaussian-mesh-momentum.py'))
    error = q(0)
    for a,b in [(0,1)]+[(2**j,2**(j+1)) for j in range(14)]:
        for j in range(256):
            l = q(Fraction(a)+Fraction((b-a)*j,256))
            r = q(Fraction(a)+Fraction((b-a)*(j+1),256))
            error += (r*r-l*l)*mesh['mesh_slope_error'](l,r)/2
    assert 8192+2*kappa<16384 and error<q('.07976')
    finite = lower-sum(net.values())*(error+2*m['tail_bound'](8192))
    print('completed matter Gaussian slope floors, continuous and finite',lower.str(12),finite.lower().str(12))
    assert finite>0
