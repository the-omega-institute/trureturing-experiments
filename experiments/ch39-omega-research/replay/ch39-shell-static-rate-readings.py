"""Exact count and minimax error for squared-amplitude data on the vertex shell."""
from math import prod
from pathlib import Path
from runpy import run_path
from flint import arb, arb_poly

root=Path(__file__).resolve().parent
model=run_path(str(root/'ch39-shell-propagator-readings.py'))
ids=[1,2,3,4,5,7,9,13,18,26,39,65,145]

def main():
    char=model['pencil'](arb(0),False).charpoly()
    den=arb_poly([v.real*(-1)**(7-i) for i,v in enumerate(char.coeffs())])
    assert all(c>0 for c in den.coeffs())
    shift=lambda i:(model['height'](i)/model['height'](1))**2
    nodes=[shift(i) for i in ids]
    variable=arb_poly([0,1])
    dual=arb_poly([0])
    amplification=arb(0)
    for i,node in enumerate(nodes):
        remaining=nodes[:i]+nodes[i+1:]
        cardinal=prod(variable-u for u in remaining)*(1/prod(node-u for u in remaining))
        coefficient=cardinal(0)*den(node)**2/den(0)**2
        assert (-1)**i*coefficient>0
        dual+=(-1)**i*cardinal*den(node)**2
        amplification+=abs(coefficient)
    assert arb('86872.34369')<amplification<arb('86872.34370')
    for i in range(1,len(model['prefix'])+1):
        if i not in ids:
            s=shift(i)
            assert abs(dual(s))<den(s)**2
    start=shift(len(model['prefix'])).lower().floor()
    for sign in (1,-1):
        tail=(den**2+sign*dual)(variable+start)
        assert all(c>0 for c in tail.coeffs())
    print('Optimal RH indices',ids,'static rate amplification',amplification.str(16))

if __name__=='__main__': main()
