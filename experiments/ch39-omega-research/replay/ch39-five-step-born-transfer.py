"""Uniform cross-type flow bound for the cubic shell's coherent five-step Born readout."""
from collections import defaultdict
from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

root = Path(__file__).parent
shell = run_path(str(root/'ch39-refined-mass-obstruction.py'))
clock = run_path(str(root/'ch39-charged-arrival-clock.py'))


def verify():
    ctx.prec = 144
    nodes, edges = shell['geometry']()
    neighbors = [[] for _ in nodes]
    for i,j,weight,charge in edges:
        assert weight>0 and abs(charge)<=1
        neighbors[i].append((j,weight))
        neighbors[j].append((i,weight))
    populations = {nodes.index(point):number for site,number in zip(clock['sites'],clock['per_site'],strict=True)
                   for point in clock['orbit'](site)}
    long = [i for i,n in populations.items() if n!=36]
    short = [i for i,n in populations.items() if n==36]
    assert sum(populations[i] for i in long)==192 and sum(populations[i] for i in short)==432
    response = {}
    for i in populations:
        row = {i:Q(1)}
        for _ in range(5):
            next_row = defaultdict(Q)
            for u,x in row.items():
                for v,y in neighbors[u]:
                    next_row[v] += x*y
            row = next_row
        response[i] = row
        assert row[i]>0
    # For g>14, each path phase has real part at least 1-25*2/785.
    cone = 1-Q(50,785)
    assert cone==Q(147,157)
    def ball(x):
        return arb(x.numerator)/x.denominator
    remaining = arb(0)
    for i in short:
        for j in long:
            # Bii^(1/2) >= cone * H(1)^5_ii at both ends.
            ratio = response[i][j]**2/(cone**2*response[i][i]*response[j][j])
            remaining += arb((populations[i]-1)*(populations[j]-1)).sqrt()*ball(ratio)
    flow = len(short)+len(long)+remaining
    assert flow < arb('56.36')
    p, q = arb(4)/13, arb(9)/13
    lag2 = p*q*(1-flow/(624*p*q))**2
    coefficient = lag2*2054195/378549
    assert lag2 > arb('.07067') and coefficient > arb('.3835')
    print('Short-to-long flow upper bound',flow.upper().str(12),'required', (624*p*q).str(12))
    print('Minimum lag-two covariance',lag2.lower().str(12),
          'coefficient error lower bound',coefficient.lower().str(12))


if __name__=='__main__':
    verify()
