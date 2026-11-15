#!/usr/bin/env python3
"""Floating point importance-sampling diagnostic, not exact evidence."""
from random import Random
from math import sqrt
from check import coefficients, terms, cas

N = 200000
rng = Random(261101492)
print('DIAGNOSTIC ONLY: floating point exponentials, seed=261101492, samples per n=', N)
print('n U_exact U_est U_se CAS_exact CAS_est CAS_se APPT_est APPT/U APPT/CAS')
for n in range(2, 9):
    d = 2*n
    c = coefficients(2,d)
    rates = [float(1-d*z/2) for z in c[:-1]]
    t = float(terms(2,d)[-1])
    hits = [0,0,0]
    for _ in range(N):
        ys = [rng.expovariate(rate) for rate in rates]
        B = -d/2*sum(float(ci)*yi for ci,yi in zip(c[:-1],ys))
        z = rng.expovariate(1)
        ed = B+z
        if ed < 0:
            continue
        x = [0.0]*d
        x[-1] = ed/d
        for j in range(d-2,-1,-1):
            x[j] = x[j+1]+ys[j]/(j+1)
        hits[0] += 1
        hits[1] += x[0] <= 3*x[-1]
        hits[2] += (x[0]-x[-2])**2 <= 4*x[-3]*x[-1]
    est = [t*h/N for h in hits]
    se = [t*sqrt((h/N)*(1-h/N)/N) for h in hits]
    u = float(sum(terms(2,d)));v = float(cas(2,d))
    # Prespecified diagnostic tolerance, not a rigorous interval.
    assert abs(est[0]-u) <= 6*se[0]
    assert abs(est[1]-v) <= 6*se[1]
    assert hits[1] <= hits[2] <= hits[0]
    print(n,*(format(a,'.9g') for a in [u,est[0],se[0],v,est[1],se[1],est[2],est[2]/u,est[2]/v]))
print('MONTE_CARLO_DIAGNOSTIC_OK')
