"""Cardinality for second-derivative matching, and a third-derivative positivity test."""
from pathlib import Path
import runpy

from flint import arb

root = Path(__file__).resolve().parent
model = runpy.run_path(str(root / 'ch39-finite-charge.py'))
q = model['q']
low = [sum(y**k*(-y).exp() for y in model['low_y']) for k in range(4)]

def bounds(boxes):
    return min(x.lower() for x in boxes).union(max(x.upper() for x in boxes))

def box(lo, hi, denominator):
    return (q(lo)/denominator).union(q(hi)/denominator)

def sample_moments():
    entropies, variances, cubic_gaps = [], [], []
    # Separate each possible left eigenvalue. The certified band width is 0.008.
    bands = [box(j, j+1, 2000) for j in range(4300, 4316)]
    high = [[y**k*(-y).exp() for k in range(4)] for y in bands]
    for j in range(5100, 5210):
        r = (-2*box(j,j+1,2000)).exp()
        deficit = [x-18*r**k*(-r).exp() for k,x in enumerate(low)]
        for i, first in enumerate(high):
            for second in high[i:]:
                w,m,b,c = [x+y+z for x,y,z in zip(deficit,first,second)]
                mean = m/w
                entropies.append(mean+w.log())
                variances.append(b/w-mean**2)
                assert m-w > 0
                cubic_gaps.append(b+(b-m)**2/(m-w)-c)
    return bounds(entropies), bounds(variances), bounds(cubic_gaps)

def binary(p):
    return -p*p.log()-(1-p)*(1-p).log()

def tail_at_entropy(eta, n):
    log_tail_count = q(n-1).log()
    # Use a floating bisection as a guess, verify both endpoints using Arb.
    import math
    left, right = 0., .5
    for _ in range(60):
        p = (left+right)/2
        trial = -p*math.log(p)-(1-p)*math.log1p(-p)+p*math.log(n-1)
        if trial < float(eta.mid()):
            left = p
        else:
            right = p
    lo, hi = q(str(left))-q('1e-13'), q(str(right))+q('1e-13')
    assert 0 < lo and hi < q('0.5')
    assert binary(lo)+lo*log_tail_count < eta < binary(hi)+hi*log_tail_count
    return lo.union(hi)

def peak_variance(eta, n):
    p = tail_at_entropy(eta,n)
    ratio = ((1-p)*(n-1)/p).log()
    return p*(1-p)*ratio**2

eta, variance, gap = sample_moments()
assert 0 < eta and eta < q(2).log()
print('entropy, target variance, minimal third moment excess')
for x in (eta,variance,gap):
    print(x.lower().str(12), x.upper().str(12))
for n in (2, 75000, 5000000000):
    lo = tail_at_entropy(eta.lower(), n)
    hi = tail_at_entropy(eta.upper(), n)
    # The derivative of the maximal variance in p has the sign of the following expression.
    assert (1-2*hi)*((1-hi)*(n-1)/hi).log() > 2
    v0, v1 = peak_variance(eta.lower(),n), peak_variance(eta.upper(),n)
    print('number of high right states',n,'attainable maximal variance',v0.lower().str(12),v1.upper().str(12))
    if n < 5000000000:
        assert v1 < variance
    else:
        assert v0 > variance
assert gap > q('4.89')
