"""Stability of the six-zero ranking under changes to the gauge matching scale."""
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

root = Path(__file__).parent
fixed = run_path(str(root/'ch39-high-zero-planck-window.py'))
q = fixed['q']


def gauge_bounds():
    centers = fixed['model']['GAUGE']
    faces = [(q(c*F(999,1000)),q(c*F(1001,1000))) for c in centers]
    box = [lo.union(hi) for lo,hi in faces]
    step = int(64*fixed['L'])
    def evaluate(initial):
        rates = []
        def observe(i, inverse, deriv):
            if i in (step-1,step):
                rates.append(deriv)
        # Endpoint enclosures at precisely L0, and rates across both neighboring slabs.
        integrate = fixed['model']['running']['integrate']
        _, start = integrate(initial,step)
        integrate(initial,step+1,observe)
        assert len(rates)==2
        return [(a,d.union(e)) for a,d,e in zip(start,*rates,strict=True)]
    complete = evaluate(box)
    sides = [[],[]]
    for i in range(3):
        for side in range(2):
            initial = box.copy()
            initial[i] = faces[i][side]
            sides[side].append(evaluate(initial)[i])
    return complete,sides


def envelopes(data, delta):
    assert abs(delta) < arb(1)/64
    complete,sides = data
    ratios = [q(x) for x in fixed['V']]
    def shift(vals):
        return [(a+delta*d)/v for (a,d),v in zip(vals,ratios,strict=True)]
    values = shift(complete)
    outer_low = max(x.lower() for x in values)
    outer_high = min(x.upper() for x in values)
    inner_low = max(x.upper() for x in shift(sides[1]))
    inner_high = min(x.lower() for x in shift(sides[0]))
    assert outer_low < inner_low < inner_high < outer_high
    coeff = (fixed['geo']['BOOK_PLANCK']/1000)**2*(3*arb(5).sqrt()/2)/(2*arb.pi())
    def height(g):
        return (coeff*(-2*(q(fixed['L'])+delta)).exp()/g-q(fixed['OFFSET'])).sqrt()
    return (height(outer_high).lower(),height(inner_high).upper(),
            height(inner_low).lower(),height(outer_low).upper())


def certified_labels(lower,upper):
    return (int(lower.zeta_nzeros().unique_fmpz())+1, int(upper.zeta_nzeros().unique_fmpz()))


def verify():
    ctx.prec = 160
    data = gauge_bounds()
    k = arb(4)/9
    h,log_two = arb(2).sqrt(),arb(2).log()
    a,b = 1+2*k*k,2*h*k
    def envelope(v):
        return (3*h-4*v)*(a+b*v).sqrt()
    def weight(g):
        return envelope((g*log_two).cos())/(g*g+arb(1)/4)
    checks = [
        (arb('-.00003').union(arb('.00003')), (11485,11490), (11486,11489),11488),
        (arb('.00009'), (11484,11488), (11484,11488),11484),
        (arb('-.0005'), (11492,11497), (11492,11496),11495)]
    error = arb('.001467')
    for delta,outer,inner,winner in checks:
        ol,il,ih,oh = envelopes(data,delta)
        assert certified_labels(ol,oh)==outer and certified_labels(il,ih)==inner
        assert inner[0]<=winner<=inner[1]
        weights = [weight(z.imag) for z in acb.zeta_zeros(outer[0],outer[1]-outer[0]+1)]
        best = weights[winner-outer[0]]
        assert all(best*(1-error)>w*(1+error) for i,w in enumerate(weights,outer[0]) if i!=winner)
        print('Offset',delta.str(8),'guaranteed winner',winner)

    # Previously certified scale-wide enclosing and sufficient bands.
    # Compare the guaranteed band champion with all members of the enclosing band.
    origin,end = 10680,11150
    values = [weight(z.imag) for z in acb.zeta_zeros(origin,end-origin+1)]
    peak = envelope((3*h*b-8*a)/(12*b))
    best = values[11107-origin]
    assert all(best>w for n,w in enumerate(values,origin) if 11077<=n and n!=11107)
    assert all(values[10693-origin]>w for n,w in enumerate(values,origin) if n!=10693)
    for rel,count,cutoff in ((arb(0),65,11129),(error,71,11148)):
        high = (peak*(1+rel)/(best*(1-rel))-arb(1)/4).sqrt()
        assert high.zeta_nzeros().unique_fmpz()==cutoff < end
        possible = []
        for n,w in enumerate(values,origin):
            if n==11107 or w*(1+rel)>best*(1-rel):
                possible.append(n)
            else:
                assert w*(1+rel)<best*(1-rel)
        assert len(possible)==count and 11488 not in possible
        print('Floating scale candidate cap:',count,'largest retained label',max(possible))


if __name__=='__main__':
    verify()
