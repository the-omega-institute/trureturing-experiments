#!/usr/bin/env python3
"""120-decimal diagnostic sweep; numerical output is not exact evidence."""
import mpmath as mp

mp.mp.dps = 120
count = samples = 0
min_real_R = mp.inf
min_outward = mp.inf
max_residual = mp.mpf(0)
max_endpoint_error = mp.mpf(0)
for n in range(8, 41, 2):
    pairs = 0
    for m in range(2, n, 2):
        if 2*(m+1) >= n:
            continue
        count += 1
        pairs += 1
        x = mp.pi/n
        theta = (m+1)*x
        h = mp.sin(m*x)/mp.sin(2*x)
        mu = mp.exp(-2j*x)
        a = sum(1-mu**s for s in range(1,m+1))/n
        b = sum(1-(-mu)**s for s in range(1,m+1))/n
        al, be = (1-mu**(m+1))/n, (1+mu**(m+1))/n
        q = 2*sum(1-mu**(2*l) for l in range(1,m//2+1))/n
        for j in range(1001):
            w = mp.mpf(j)/1000
            R = 1+w/h*mp.exp(-1j*(m+2)*x)+w*w/(4*h*h)*mp.exp(2j*m*x)
            ell = (a+b+w/n-(b-a)*mp.sqrt(R))/2
            min_real_R = min(min_real_R,mp.re(R))
            assert mp.re(R) >= mp.mpf(3)/4
            residual = (ell-a)*(ell-b)-w*(be*(ell-a)+al*(ell-b))/2
            max_residual = max(max_residual,abs(residual))
            assert abs(residual) < mp.mpf('1e-115')
            if j == 0:
                assert abs(ell-a) < mp.mpf('1e-115')
            elif j == 1000:
                max_endpoint_error = max(max_endpoint_error,abs(ell-q))
                assert abs(ell-q) < mp.mpf('1e-115')
            else:
                distance = -mp.re(mp.exp(1j*theta)*(ell-a))
                assert distance > mp.mpf('1e-100')
                min_outward = min(min_outward,distance)
                # Check every vertex lies on the inward side of the chosen edge.
                # A convex-hull point then cannot equal this exterior eigenvalue.
                samples += 1
        for k in range(n):
            vertex = sum(1-mu**s for s in range(1,k+1))/n
            for v in [vertex,mp.conj(vertex)]:
                assert -mp.re(mp.exp(1j*theta)*(v-a)) <= mp.mpf('1e-115')
    print(f"n={n}: {pairs} admissible m, 1001 grid points each OK")
print("NUMERICAL_DIAGNOSTIC: dps=120; grid w=j/1000, j=0,...,1000")
print("pairs",count,"strict-interior samples",samples)
print("min Re R",mp.nstr(min_real_R,35))
print("min positive outward distance",mp.nstr(min_outward,35))
print("max characteristic residual",mp.nstr(max_residual,10))
print("max endpoint error",mp.nstr(max_endpoint_error,10))
print("HIGH_PRECISION_CHECKS_OK")
