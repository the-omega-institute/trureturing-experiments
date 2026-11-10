#!/usr/bin/env python3
"""Exact finite checks for the Gutkovskiy nested floor recurrences.

This program checks finite data and algebra, not convergence at infinity.
Only integer arithmetic, Fraction, and SymPy exact algebra enter assertions.
"""
from fractions import Fraction
from hashlib import sha256
from math import isqrt
import sys
import sympy as sp

assert sys.version_info >= (3, 12)
N = 1_000_000
SOURCE = {
    3: ("A366860", "1640df58179f288b880ac85d99179c6522db5ecdfe34feb14bc65fe2b12cd348",
        [0,1,2,3,3,4,5,6,6,7,8,9,9,10,11,12,13,14,15,15,16,17,18,18,19,20,21,21,22,23]),
    4: ("A366870", "611ec1e8339fc5db5e2402a11c7019c39dd44a3e884b9d9ae464933f9f318934",
        [0,1,2,3,4,4,5,6,7,8,8,9,10,11,12,12,13,14,15,16,16,17,18,19,20,21,22,23,24,24]),
    5: ("A366871", "9865622a033400a4ce893d0d792c1c85b5fa6ab1719c6d0c8aa7c365af807e92",
        [0,1,2,3,4,5,5,6,7,8,9,10,10,11,12,13,14,15,15,16,17,18,19,20,20,21,22,23,24,25]),
}


def generate(m, limit):
    a = [0] * (limit + 1)
    for n in range(1, limit + 1):
        k = a[n - 1] // m
        assert 0 <= k < n
        a[n] = n - a[k]
        assert 0 <= a[k] <= k
        assert 0 <= a[n] <= n
        assert m * a[n] >= (m - 1) * n + 1
        assert m * m * k >= (m - 1) * (n - 1) - m * m
        if n >= 2:
            # |k/n - (a[n-1]/(n-1))/m| <= (1+1/m)/n.
            assert abs(m * (n - 1) * k - n * a[n - 1]) <= (m + 1) * (n - 1)
        assert a[n] + a[k] == n
    return a


def root_bracket(m):
    scale = 10**18
    d = m*m + 4*m
    q = isqrt(d * scale**2)
    assert q*q <= d*scale**2 < (q+1)**2
    lo, hi = Fraction(q-m*scale, 2*scale), Fraction(q+1-m*scale, 2*scale)
    assert lo > 0
    assert lo*lo + m*lo - m <= 0 < hi*hi + m*hi - m
    return lo, hi


print("Exact arithmetic: Python >= 3.12; SymPy", sp.__version__)
print("Source: Gutkovskiy OEIS NAME/FORMULA; Eldar b-files, n=0..10000.")
print("Source digests hash comma-separated integer values followed by one newline.")
closed = {3: (sp.sqrt(21)-3)/2, 4: 2*(sp.sqrt(2)-1), 5: (3*sp.sqrt(5)-5)/2}
for m, c in closed.items():
    uniform = (sp.sqrt(m*m+4*m)-m)/2
    assert sp.simplify(c-uniform) == 0
    assert sp.simplify(c*c+m*c-m) == 0
    assert c.is_positive is True
    # Two-digit period [0; 1,m,1,m,...] has this positive fixed point.
    assert sp.simplify(c - 1/(1+1/(m+c))) == 0
    lower, upper = Fraction(0), Fraction(1)
    lo, hi = root_bracket(m)
    for _ in range(30):
        def t(x): return 1/(1+1/(m+x))
        lower, upper = t(lower), t(upper)
        assert lower < upper
    assert lower*lower + m*lower - m < 0 < upper*upper + m*upper - m
    assert lower <= hi and lo <= upper
    assert upper - lower < Fraction(1, 10**25)
    print(f"m={m}: closed form={c}; quadratic and periodic-CF identities exact")
    print(f"  root bracket=[{lo},{hi}]; 30-period rational CF bracket=[{lower},{upper}]")
    a = generate(m, N)
    anum, digest, prefix = SOURCE[m]
    assert a[:len(prefix)] == prefix
    computed = sha256((",".join(map(str, a[:10001]))+"\n").encode()).hexdigest()
    assert computed == digest
    print(f"  FIDELITY {anum}: all 10001 b-file values agree; values_sha256={computed}")
    for n in (10,100,1000,10000,100000,N):
        print(f"  n={n}: a(n)={a[n]}; a(n)/n={Fraction(a[n],n)}")
    # A finite tail statement only; root brackets make the error bound exact.
    n_min = n_max = 10000
    for n in range(10001, N+1):
        if a[n]*n_min < a[n_min]*n: n_min = n
        if a[n]*n_max > a[n_max]*n: n_max = n
    rmin, rmax = Fraction(a[n_min],n_min), Fraction(a[n_max],n_max)
    error = max(hi-rmin,rmax-lo)
    assert error < Fraction(1,1000)
    residual_bound = Fraction(abs(a[N]**2+m*a[N]*N-m*N*N),m*N*N)
    assert max(hi-Fraction(a[N],N),Fraction(a[N],N)-lo) <= residual_bound
    print(f"  finite tail 10000..{N}: min={rmin} at n={n_min}; max={rmax} at n={n_max}")
    print(f"  |a(n)/n-c| < 1/1000 throughout that finite tail; certified upper bound={error}")
    print(f"  endpoint root-error upper bound from quadratic residual={residual_bound}")

for m in range(3,41):
    generate(m,5000)
print("General-family structural checks: m=3..40, n=0..5000.")
x,y,m = sp.symbols('x y m')
assert sp.expand((y*y-x*x)-(y-x)*(y+x)) == 0
assert sp.expand((1-x*x/m)-(1-y*y/m)-(y*y-x*x)/m) == 0
assert sp.expand((x*x+m*x-m)-(y*y+m*y-m)-(x-y)*(x+y+m)) == 0
print("Oscillation subtraction, factorization, and root-difference identities: exact.")
print("Finite checks do not prove convergence; the all-index argument remains for the later probe.")
print("ALL_CHECKS_OK")
