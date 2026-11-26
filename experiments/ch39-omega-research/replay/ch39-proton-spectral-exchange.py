"""Sharp RH phase-cancellation bounds for a specified tree-level triplet exchange."""
from fractions import Fraction as Q
from math import prod
from pathlib import Path
from runpy import run_path

root = Path(__file__).resolve().parent
old = run_path(str(root/'ch39-eight-time-pole-control.py'))
B, node, multiply = (old[k] for k in ('B','node','mul'))
sqrt_bounds = old['root_cell']
time_error = Q(1, 2**272)


def legs(table, x):
    nodes, weights = table
    return [w*z/(B(1)+B(x)*z) for z,w in zip(nodes,weights,strict=True)]


def floor_bounds(table, x):
    values = legs(table,x)
    tail = sum(values[1:],B(0))
    gap = values[0]-tail
    return gap, gap/(values[0]+tail)


def check_extrema():
    first, second = node(1), node(2)
    mass_floor = old['base']['data']('ch39-natural-mass-cell.json')
    mass = B(Q(mass_floor,2**288),Q(mass_floor+1,2**288))
    # Equal exchange lengths on the first two heights, with opposite signs.
    low_first = second/(first+second)
    assert (first/mass).hi < low_first.lo
    pairs = [1,2,3,4,8,9,51,52]
    assert all(i+1==j for i,j in zip(pairs[::2],pairs[1::2]))
    for ids,raw,relative in (([1,2,0],'.00028112738223','.17473815523'),
                              (pairs,'.00055971852907','.34789853304')):
        table = old['quadrature'](ids,time_error)
        nodes,weights = table
        squares = tuple(j*j for j in range(1,len(ids)))
        at_first = multiply(*(B(1)+B(a)*first for a in squares))
        derivative = multiply(*(first-v for v in nodes[1:]))
        slope = at_first/(first*derivative)
        # In both rules H/W is nonnegative on all unlisted genuine heights.
        # On 0 <= x <= max(squares), the unsigned positive error coefficient
        # in the total length is at most max(squares)**len(squares).
        last = squares[-1]
        if len(squares)==7:
            assert (first/(B(1)+B(last)*first)*slope).lo > last**len(squares)
        else:
            assert (first*slope).lo > prod(squares)
        gap,ratio = floor_bounds(table,0)
        assert Q(raw) < gap.lo <= gap.hi < Q(raw)+Q(1,10**14)
        assert Q(relative) < ratio.lo <= ratio.hi < Q(relative)+Q(1,10**11)
        print(len(squares),'nonconstant readings, raw floor',float(gap.lo),'relative floor',float(ratio.lo))
    table = old['quadrature'](pairs,time_error)
    cumulative = table[1][0]
    for weight in table[1][1:]:
        assert cumulative.lo>weight.hi
        cumulative = cumulative+weight
    x_lo, x_hi = Q('496.77998547'), Q('496.77998548')
    assert floor_bounds(table,x_lo)[0].lo>0>floor_bounds(table,x_hi)[0].hi
    assert 22**2 < x_lo < x_hi < 23**2
    print('The first spacelike cancellation on the exact eight-row problem is bracketed by',float(x_lo), float(x_hi))


def check_benchmark():
    # The mass and the seven-site exact table come from the same threshold fit.
    table = old['quadrature']([1,2,3,4,11,12,3200],Q(0))
    gap,_ = floor_bounds(table,0)
    mf = B(69169689551423,69169689551424)
    q02 = old['base']['physical_inputs']()[2]
    assert (mf**2/table[0][0]).lo > q02.hi
    # A rational enclosure for pi, or any tighter enclosure, suffices here.
    pi = B(Q(333,106),Q(355,113))
    second_per_year = 31557600
    width = Q('6.582119569e-25')/(Q('2.4e34')*second_per_year)
    mp, mpi = Q('.938272'), Q('.134977')
    amplitude = sqrt_bounds(B(32*width/(mp*(1-(mpi/mp)**2)**2))*pi)
    upper = amplitude*mf**2/gap
    assert Q('.08422') < upper.lo <= upper.hi < Q('.08424')
    # With |h_D h_L W_eff|=.10 GeV^2, even the greatest lifetime falls short.
    limit_fraction = (upper/B(Q('.1')))**2
    assert limit_fraction.hi < Q('.710')
    print('Benchmark upper allowed |h_D h_L W_eff| / GeV^2:', float(upper.hi))


if __name__ == '__main__':
    check_extrema()
    check_benchmark()
