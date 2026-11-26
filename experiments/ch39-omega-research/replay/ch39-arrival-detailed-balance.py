"""Detailed-balance test for the arithmetic four-tap comparison."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
from runpy import run_path


def coefficient_identity():
    p, q = F(4,13), F(9,13)
    gains1, gains2 = [-1/q,1/p], [-2/q,1/p]
    target = [F(x,13) for x in (18,-19,-3,4)]
    v = [F(x,2054195) for x in (11259,64962,-177696,124632)]
    assert sum(a*b for a,b in zip(v,target,strict=True))==0
    intercept, marginal = F(-278,12155), F(9287,24310)
    phi = {(0,0):F(-32,65),(0,1):F(-274,715),(1,0):F(-61,110),(1,1):0}
    for i,j,k in product(range(2),repeat=3):
        assert i*k == (v[i+j+k]*gains1[i]*gains2[j]+intercept+marginal*i
                       +phi[i,j]-phi[j,k])
    assert intercept + marginal*p == p**2
    assert (p*q)/sum(abs(x) for x in v)==F(48620,42061)


def symmetric_kernel():
    clock = run_path(str(Path(__file__).with_name('ch39-charged-arrival-clock.py')))
    offsets = clock['routes']()
    # Include the reverse path, with charge reversed, for every saved offset.
    shifts = {(i,j):{a,(-offsets[j,i])%3} for (i,j),a in offsets.items()}
    assert all({(-o)%3 for o in row}==shifts[j,i] for (i,j),row in shifts.items())
    sizes = clock['per_site']
    for i in range(6):
        to_type = [F(0),F(0)]
        for j in range(6):
            if i<4 and j<4:
                factor = 2 if (i<2)!=(j<2) else 0
            else:
                factor = 1
            if factor:
                assert (i,j) in shifts
                # Sum over target channels and allowed coordinate directions.
                per_channel = F(factor,208*len(shifts[i,j]))
                to_type[int(j>=4)] += sizes[j]*len(shifts[i,j])*per_channel
        assert to_type==[F(9,13),F(4,13)]
    print('Symmetric five-step transition; from every channel type probabilities are 9/13 and 4/13.')


if __name__=='__main__':
    coefficient_identity()
    symmetric_kernel()
