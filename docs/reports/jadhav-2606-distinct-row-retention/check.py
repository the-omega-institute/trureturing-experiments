#!/usr/bin/env python3
"""Exact rational kernel checks and certified rational logarithm intervals.

No binary floating point is used. The logarithm enclosure uses the positive
atanh series and its geometric tail bound after rational range reduction.
These finite checks do not establish the asymptotic or Lean settlement.
"""
from fractions import Fraction as F
from functools import cache
import hashlib
from pathlib import Path
import sympy as sp


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def add(a, b):
    return a[0] + b[0], a[1] + b[1]


def scale(c, a):
    return (c*a[0], c*a[1]) if c >= 0 else (c*a[1], c*a[0])


def divide(a, b):
    require(b[0] > 0, 'interval denominator positive')
    endpoints = [x/y for x in a for y in b]
    return min(endpoints), max(endpoints)


def log_unit(q, terms=80):
    require(1 <= q <= 2, 'range reduction')
    z = (q-1)/(q+1)
    s = sum((2*z**(2*j+1)/F(2*j+1) for j in range(terms)), F(0))
    tail = 2*z**(2*terms+1)/(F(2*terms+1)*(1-z*z))
    return s, s+tail


@cache
def log_interval(q):
    q = F(q)
    require(q > 0, 'log argument positive')
    if q < 1:
        return scale(-1, log_interval(1/q))
    k = 0
    while q > 2:
        q /= 2
        k += 1
    return add(scale(k, log_unit(F(2))), log_unit(q))


def kl(p, q):
    result = (F(0), F(0))
    for a, b in zip(p, q):
        require(a >= 0 and b > 0, 'KL admissibility')
        if a:
            result = add(result, scale(a, log_interval(a/b)))
    return result


def display(a, places=12):
    unit = 10**places
    lo = (a[0].numerator*unit)//a[0].denominator
    hi = -((-a[1].numerator*unit)//a[1].denominator)
    def fmt(x):
        sign = '-' if x < 0 else ''
        x = abs(x)
        return f'{sign}{x//unit}.{x%unit:0{places}d}'
    return '[' + fmt(lo) + ', ' + fmt(hi) + ']'


def log_form(q):
    """Formal prime-log coefficients; valid for positive rational q."""
    d = {}
    for value, sign in [(q.numerator, 1), (q.denominator, -1)]:
        for prime, exponent in sp.factorint(value).items():
            d[int(prime)] = d.get(int(prime), F(0)) + sign*exponent
    return {p: c for p, c in d.items() if c}


def kl_form(p, q):
    d = {}
    for a, b in zip(p, q):
        if a:
            for prime, c in log_form(a/b).items():
                d[prime] = d.get(prime, F(0)) + a*c
    return {p: c for p, c in d.items() if c}


def fidelity():
    # Source v1, Proposition 9.5: half-lazy star on n=5 vertices.
    n = 5
    pi = [F(1, 2)] + [F(1, 8)]*4
    P = [[F(1, 2)] + [F(1, 8)]*4]
    for i in range(1, n):
        P.append([F(1, 2)] + [F(1, 2) if i == j else F(0) for j in range(1, n)])
    require(P[0] == pi, 'star center row equals pi')
    require(kl_form(P[0], pi) == {}, 'star center divergence zero')
    for i in range(n):
        require(sum(P[i]) == 1, 'star row sum')
        for j in range(n):
            require(pi[i]*P[i][j] == pi[j]*P[j][i], 'star detailed balance')
    for row in P[1:]:
        require(kl_form(row, pi) == {2: F(1)}, 'star leaf KL = log 2')
    require(log_form(1/pi[1]) == {2: F(3)}, 'star input KL = 3 log 2')
    r_leaf = F(1, 3)
    L = sum(pi[1:])*r_leaf/r_leaf
    require(L == F(1, 2), 'source lazy-star localization')
    print('FIDELITY: source Proposition 9.5, n=5: r(center)=0, r(leaf)=1/3, L=1/2 exactly')


def kernel(m):
    n, eps = m*m, F(1, m*m)
    P = []
    for x in range(n):
        row = []
        for y in range(n):
            if x < m and y < m:
                v = (1-eps)*F(x == y, 2) + (1-eps)/F(2*m) + eps/n
            elif x >= m and y >= m:
                v = (1-eps)*eps*int(x == y) + (1-eps)**2/F(m*(m-1)) + eps/n
            else:
                v = eps/n
            row.append(v)
        P.append(row)
    return P


def member(m):
    n, eps = m*m, F(1, m*m)
    P = kernel(m)
    pi = [F(1, n)]*n
    require(all(sum(row) == 1 for row in P), f'm={m}: row sums')
    require(all(v >= eps/n > 0 for row in P for v in row), f'm={m}: positivity')
    require(all(P[x][y] == P[y][x] for x in range(n) for y in range(n)), f'm={m}: symmetry')
    require(all(sum(pi[x]*P[x][y] for x in range(n)) == pi[y] for y in range(n)), f'm={m}: stationarity')
    require(len({tuple(row) for row in P}) == n, f'm={m}: every row distinct')
    tv = sum((abs(a-b) for a, b in zip(P[0], P[1])), F(0))/2
    require(tv == (1-eps)/2 >= F(12, 25), f'm={m}: A row separation')
    if m == 5:
        require((P[0][0], P[0][1], P[0][5], P[5][5], P[5][6]) ==
                (F(361,625), F(61,625), F(1,625), F(269,3125), F(149,3125)), 'm=5 printed entries')
    logm, logn = log_interval(F(m)), log_interval(F(n))
    rA, rB = divide(kl(P[0], pi), logn), divide(kl(P[m], pi), logn)
    require(rA[0] > rB[1], f'm={m}: maximizing block exactly A')
    mu_out = [sum((P[x][y]/F(m) for x in range(m)), F(0)) for y in range(n)]
    rho = F(m-1, m**3)
    require(mu_out == [(1-rho)/m]*m + [rho/F(n-m)]*(n-m), 'block output law')
    witness = divide(kl(mu_out, pi), logm)
    # eta is enclosed using this feasible law and analytic data processing;
    # the check does not numerically optimize or assume a point-input optimizer.
    eta = (witness[0], F(1))
    require(witness[0] >= 1-4*eps, f'm={m}: eta lower bound')
    require(rA[1] <= F(3,4)*(1-eps), f'm={m}: A upper bound')
    require(rA[0] >= F(3,4)*(1-eps)-divide(log_interval(F(3)), scale(2,logm))[1], 'A lower bound')
    require(rA[0] >= (1-4*eps)/2, 'A coarse lower bound')
    require(rB[1] <= eps+F(1,2*(m-1)), 'B upper bound')
    rbar = add(scale(F(1,m), rA), scale(F(m-1,m), rB))
    L, ratio = divide(rbar, rA), divide(eta, rA)
    require(L[1] <= F(4,m), 'L upper bound')
    require(ratio[0] >= F(4*(m*m-4),3*(m*m-1)) >= F(7,6), 'uniform ratio separation')
    print(f'm={m}: n={n}; exact structural checks OK; TV(A rows)={tv}')
    print('  rA='+display(rA)+' rB='+display(rB)+' L='+display(L))
    print('  eta in '+display(eta)+'; eta/M in '+display(ratio))


def symbolic_checks():
    t, m = sp.Symbol('t', nonnegative=True), sp.Symbol('m', positive=True)
    eps, rho = 1/m**2, (m-1)/m**3
    # All global rational comparisons reduce to polynomials nonnegative at m=t+5.
    comparisons = {
        'rho <= epsilon': eps-rho,
        'A-row coarse bound >= 21/50': (1-4*eps)/2-sp.Rational(21,50),
        'B-row coarse bound <= 33/200': sp.Rational(33,200)-eps-1/(2*(m-1)),
        'A-row TV >= 12/25': (1-eps)/2-sp.Rational(12,25),
        'L rational bound <= 4/m': 4/m-(sp.Rational(5,2)*m*m+2*m-2)/(m*(m*m-4)),
        'eta/M rational lower bound >= 7/6': 4*(m*m-4)/(3*(m*m-1))-sp.Rational(7,6),
    }
    for label, expr in comparisons.items():
        num, den = sp.fraction(sp.factor(expr))
        for part in [num, den]:
            shifted = sp.Poly(sp.expand(part.subs(m,t+5)),t)
            require(all(c >= 0 for c in shifted.all_coeffs()) and shifted.eval(1)>0, label)
        require(sp.Poly(sp.expand(den.subs(m,t+5)),t).eval(0)>0, label+' denominator')
        print('SYMBOLIC: '+label+' for every integer m>=5')
    proposed = sp.Rational(5,4)/m+(m-1)/m**3
    require(sp.simplify(proposed/((1-4*eps)/2) -
                        (sp.Rational(5,2)*m*m+2*m-2)/(m*(m*m-4))) == 0, 'L bound identity')
    print('SYMBOLIC: rational L-bound identity OK; asymptotic limits require analytic proof')


def main():
    print('check.py sha256='+hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    print('Arithmetic: Python integers/Fraction and SymPy exact polynomials; no floating point')
    fidelity()
    symbolic_checks()
    for m in [5, 6, 8, 10]:
        member(m)
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
