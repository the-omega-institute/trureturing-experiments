#!/usr/bin/env python3
"""Herrera, arXiv:2311.17448v1, Conjecture 2 (Conjecture 1.10 in print): for every non-negative operator monotone f on
[0,inf) and every c > 0, E_f(c) <= f(c), where
E_f(c) = inf over g' in M_+ (Fourier transforms of finite positive Borel measures on R) of
         sup_{x>=0}(f-g)(x) - inf_{x>=0}(f-g)(x) + c g'(0).

Refutation route: for every admissible g, (N^3+N+2) e >= D_N(f,c) := (N+1) f(2c) - (N+2) f(c) - f(Nc), where e = J - f(c).
Checks:
 (1) exact (sympy): D_N(f_a,c) = r^2 (N^2-(4r+3)N-2) / ((1+r)(1+2r)(1+Nr)) for f_a(x) = x/(a+x), r = c/a;
 (2) exact: N = 11, c = a gives D = 7/12 and the gap 1/2304, i.e. E_{f_a}(a) >= 1153/2304 > 1/2 = f_a(a);
 (3) N >= 4r+4 gives a positive numerator, for r on a grid;
 (4) numeric: 1 - cos(N u) <= N^2 (1 - cos u) on a grid, N = 1..12;
 (5) numeric consistency (numpy only): the source's Gaussian family g'(x) = A exp(-b x^2) (the Fourier transform of a positive
     Gaussian measure) gives J >= 1153/2304 on a parameter grid for f_1 at c = 1, and its best J is near the source's
     constant C < 1.01975 times f(1) = 1/2."""
import sys, math
import sympy as sp
import numpy as np
erf = np.vectorize(math.erf)
bad = 0
def report(label, ok):
    global bad
    print(label, ok); bad += 0 if ok else 1
N, r, a, c = sp.symbols('N r a c', positive=True)
f = lambda x: x / (a + x)
D = (N + 1) * f(2 * c) - (N + 2) * f(c) - f(N * c)
target = r**2 * (N**2 - (4 * r + 3) * N - 2) / ((1 + r) * (1 + 2 * r) * (1 + N * r))
report('(1) D_N(f_a,c) closed form', sp.simplify(D.subs(c, r * a) - target) == 0)
D11 = sp.nsimplify(D.subs({N: 11, c: a}))
gap = D11 / (11**3 + 11 + 2)
report(f'(2) N=11, c=a: D={D11}, gap={gap}, bound={sp.Rational(1,2)+gap}', D11 == sp.Rational(7, 12) and gap == sp.Rational(1, 2304))
ok3 = True
for rv in [x / 8 for x in range(1, 160)]:
    Nv = math.ceil(4 * rv + 4)
    ok3 &= Nv**2 - (4 * rv + 3) * Nv - 2 >= Nv - 2 > 0
report('(3) positive numerator for N >= 4r+4', ok3)
u = np.linspace(-20, 20, 200001)
ok4 = all(np.all(1 - np.cos(n * u) <= n**2 * (1 - np.cos(u)) + 1e-12) for n in range(1, 13))
report('(4) 1-cos(Nu) <= N^2 (1-cos u)', ok4)
x = np.linspace(0, 400, 400001)
fx = x / (1 + x)
best = 10.0
for A in np.linspace(0.05, 1.0, 96):
    for b in np.geomspace(0.01, 10, 120):
        g = A * math.sqrt(math.pi / (4 * b)) * erf(math.sqrt(b) * x)
        h = fx - g
        J = h.max() - h.min() + A
        best = min(best, J)
report(f'(5) Gaussian family: best J = {best:.6f} >= 1153/2304 = {1153/2304:.6f} and <= 1.01975/2', 1153 / 2304 <= best <= 1.01975 / 2 + 1e-3)
print('bad=', bad)
sys.exit(1 if bad else 0)
