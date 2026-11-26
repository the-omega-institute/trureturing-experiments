"""A finite eight-site counterexample for offsets 1,...,6 and every J >= 232."""
from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path

root = Path(__file__).resolve().parent
old = run_path(str(root/'ch39-eight-time-pole-control.py'))
B, node, mul, sqrt_cell, base = (old[k] for k in ('B','node','mul','root_cell','base'))
common = [1, 2, 3, 4, 11, 12]
limit, nearby = [old['quadrature'](common+[i], Q(0)) for i in (0, 115)]
scale = B(Q('.000835'))


def check_far_readings():
    # Enclose squared heights on a binary grid, then sum outward-rounded
    # reciprocal terms. All higher nodes contribute the certified tail mass.
    height_grid, sum_grid = 2**40, 2**72
    heights = [(h*h*height_grid//2**800, (h+1)**2*height_grid//2**800+1)
               for h, *_ in old['prefix']]
    def partial(j):
        shift = j*j*height_grid
        numerator = sum_grid*height_grid
        lo = sum(numerator//(h+shift) for _, h in heights)
        hi = sum(numerator//(l+shift)+1 for l, _ in heights)
        return B(Q(lo,sum_grid), Q(hi,sum_grid))
    cell = base['data']('ch39-natural-mass-cell.json')
    mass = B(Q(cell, 2**288), Q(cell+1, 2**288))
    tail = B(1)-partial(0)/mass
    assert tail.lo > 0
    def natural(j):
        return partial(j)/mass + tail/B(1, 1+j*j*node(len(heights)).hi)
    def reading(table,j):
        return sum((w/(B(1)+B(j*j)*v) for v,w in zip(*table)), B(0))
    assert natural(231).hi < reading(nearby,231).lo
    assert natural(232).hi < Q('.17')
    for j in range(232, 480):
        assert natural(j).lo > reading(nearby,j).hi
    moment = sum((w/v for v,w in zip(*nearby)),B(0))
    assert 480**2*natural(480).lo > moment.hi
    # The mass on the free node increases as its coordinate rises from 0.
    # A genuine height at least 2J keeps at least 4/5 of that contribution.
    assert Q(4,5)*limit[1][-1].lo > Q('.17')
    print('The two seven-site tables straddle the natural seventh reading for every J >= 232.')


def check_preparations():
    z, low = limit
    high = nearby[1]
    fixed = [B(min(a.lo,b.lo), max(a.hi,b.hi)) for a,b in zip(low[:6],high[:6])]
    edges = [sqrt_cell(w*v) for w,v in zip(fixed,z[:6])]
    first = edges[0]
    other = max([x.hi for x in edges[1:]]+[sqrt_cell(node(115)).hi])
    length = sum(v.lo for v in edges)
    assert first.lo > other and first.lo**2 > scale.hi
    assert first.hi**2*old['b41'].hi < scale.lo
    cosines = old['cosines']
    loss = 1-sqrt_cell((B(1)+cosines[41])/B(2)).lo
    assert length > 2*first.hi+other*loss and length > 2*other+first.hi*loss
    sine = [B(0)]+[sqrt_cell((B(1)-v)/B(2)) for v in cosines[1:]]
    assert all(first.lo*(sine[i]-sine[i-1]).lo > other*sine[1].hi for i in range(1,42))


def check_crossing():
    # Monotonicity of the interpolating rational remainder in the moving
    # node fixes the larger response in each pole gap.
    worst = (nearby,limit,nearby)
    assert old['boundary'](scale,worst,False).hi < Q('.262')
    p, modes, q2, target = base['physical_inputs']()
    high = B(650)
    poles = [B(1)/node(i) for i in range(1,5)]
    assert poles[2].hi < 650 < poles[3].lo and (scale*high).hi<8
    H = B(125**2)*high/q2
    assert ((B(1)-p)/(B(1)-H*node(1))).hi < 1
    def secular(y,mode,table):
        r = sum((w/(B(1)-B(y)*v) for v,w in zip(*table)),B(0))
        return scale*(B(y)-H)-mode+mode*(B(1)-p)*r
    upper = [[2060,1980,1332,2036,1590,2016],
             [339,336,299,338,315,337], [569,566,526,568,546,568]]
    exposure = mul(*((high/t)**25 for t in poles[:3]))
    for j,(row,small) in enumerate(zip(upper,(limit,nearby,limit),strict=True)):
        for raw,mode,n in zip(row,modes,base['MODES'],strict=True):
            y = Q(raw,10**6 if j==0 else 1)
            assert (poles[j-1].hi if j else 0)<y<poles[j].lo
            assert secular(y,mode,small).lo>0
            exposure = mul(exposure,(high/B(y))**n)
    assert (exposure**11/target).lo > 1
    for mode in modes:
        assert secular(Q(650),mode,limit).hi < 0
        y = Q(1,10**5)
        assert (scale*B(y)-mode+mode*(B(1)-p)/(B(1)-B(y)*node(1))).hi < 0
        assert (q2/high).lo*y > 1000**2
    print('For any such seventh offset, the crossing has 144 open doublets at kinetic scale 0.000835.')


if __name__ == '__main__':
    check_far_readings()
    check_preparations()
    check_crossing()
