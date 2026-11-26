"""Localize the floating gauge matching interval using endpoint and derivative enclosures."""
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

base = run_path(str(Path(__file__).with_name('ch39-prime-window-matching-scale.py')))
fixed, q = base['fixed'], base['q']


def trajectory(initial, check_outside, tolerance=0):
    steps = 64
    center = int(steps*fixed['L'])
    radius = steps//8
    derivatives = None
    ratios = [q(v) for v in fixed['V']]
    def observe(i, inverse, rate):
        nonlocal derivatives
        if center-radius <= i < center+radius:
            derivatives = rate if derivatives is None else [a.union(b) for a,b in zip(derivatives,rate)]
        elif check_outside and i >= 23*steps:
            normalized = [x/v for x,v in zip(inverse,ratios)]
            assert max(x.lower() for x in normalized)>min(x.upper() for x in normalized)+2*tolerance
    integrate = fixed['model']['running']['integrate']
    _, at_center = integrate(initial,center)
    integrate(initial,33*steps if check_outside else center+radius,observe)
    return list(zip(at_center,derivatives,strict=True))


def scale_bounds(tolerance=0):
    faces = [(q(c*F(999,1000)),q(c*F(1001,1000))) for c in fixed['model']['GAUGE']]
    box = [a.union(b) for a,b in faces]
    all_data = trajectory(box,True,tolerance)
    side_data = [[],[]]
    for i in range(3):
        for side in range(2):
            initial = box.copy()
            initial[i] = faces[i][side]
            side_data[side].append(trajectory(initial,False)[i])
    return all_data, side_data


def bands(data, tolerance=0):
    scale = 65536
    ratios = [q(x) for x in fixed['V']]
    coeff = (fixed['geo']['BOOK_PLANCK']/1000)**2*(3*arb(5).sqrt()/2)/(2*arb.pi())
    constant = coeff*(-2*q(fixed['L'])).exp()
    offset = q(fixed['OFFSET'])
    def shifted(row, delta):
        return [(a+delta*b)/v for (a,b),v in zip(row,ratios,strict=True)]
    outer,inner = [],[]
    for i in range(-scale//8,scale//8):
        left, right = arb(i)/scale,arb(i+1)/scale
        slab = left.union(right)
        all_values = shifted(data[0],slab)
        low = (max(v.lower() for v in all_values)-tolerance).lower()
        high = (min(v.upper() for v in all_values)+tolerance).upper()
        def height(dt,g):
            return (constant*(-2*dt).exp()/g-offset).sqrt()
        if low > high:
            continue
        outer.append((height(right,high).lower(),height(left,low).upper()))
        low = max(v.upper() for v in shifted(data[1][1],slab))
        high = min(v.lower() for v in shifted(data[1][0],slab))
        if low < high:
            inner.append((height(right,high).upper(),height(left,low).lower()))
    inner.sort(key=lambda pair:float(pair[0]))
    lower, covered = inner[0]
    for lo,hi in inner[1:]:
        assert lo < covered
        covered = max(hi,covered)
    return (min(x for x,_ in outer),max(y for _,y in outer)), (lower,covered)


def verify():
    ctx.prec = 160
    outer,inner = bands(scale_bounds())
    labels = [base['certified_labels'](*ends) for ends in (outer,inner)]
    assert labels==[(10815,12092),(10905,12054)]
    print('outer, inner',labels,[[float(x) for x in ends] for ends in (outer,inner)])
    first, guaranteed = labels[0][0], labels[1][0]
    h, k, log_two = arb(2).sqrt(),arb(4)/9,arb(2).log()
    a,b = 1+2*k*k,2*h*k
    W = lambda v: (3*h-4*v)*(a+b*v).sqrt()
    peak = W((3*h*b-8*a)/(12*b))
    def weight(g):
        return W((g*log_two).cos())/(g*g+arb(1)/4)
    initial = [weight(v.imag) for v in acb.zeta_zeros(guaranteed,25)]
    bound = max(initial)
    eps = arb('.001467')
    last = int((peak*(1+eps)/(bound*(1-eps))-arb(1)/4).sqrt().zeta_nzeros().unique_fmpz())
    values = [weight(v.imag) for v in acb.zeta_zeros(first,last-first+1)]
    win = max(range(guaranteed-first,len(values)),key=lambda j:float(values[j]))
    assert first+win==10908 and last<labels[1][1]
    for e in (arb(0),eps):
        kept = [first+j for j,w in enumerate(values) if j==win or w*(1+e)>values[win]*(1-e)]
        expected = [10843,10854,10871,10875,10908]+([10918] if e>0 else [])
        assert kept==expected
        assert all(w*(1+e)<values[win]*(1-e) for j,w in enumerate(values) if first+j not in kept)
        print('guaranteed best, contenders',first+win,kept)


if __name__=='__main__':
    verify()
