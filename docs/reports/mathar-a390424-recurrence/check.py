#!/usr/bin/env python3
"""Exact checks for OEIS A390424 and Mathar's eventual recurrence.

Requires Python 3 and SymPy. No floating point arithmetic is used.
The sequence is constructed from g=1+x*g^3, independently of the recurrence.
"""
from fractions import Fraction
from math import comb
import sympy as s

N = 80
SOURCE = [
    1, 4, 26, 176, 1209, 8360, 57998, 403048, 2803322, 19505448,
    135732626, 944461632, 6570687403, 45702062144, 317791588456,
    2209127828768, 15352038918434, 106653736269992, 740714729403562,
    5142712676327408, 35694581962366591, 247676081665181864,
    1718070166963868430,
]

def mul(a, b):
    return [sum(a[j] * b[i-j] for j in range(i+1)) for i in range(N+1)]

# Each g[n] depends only on g[0],...,g[n-1].
g = [1] + [0] * N
for n in range(1, N+1):
    g[n] = sum(g[i] * g[j] * g[n-1-i-j]
               for i in range(n) for j in range(n-i))
assert g == [Fraction(comb(3*n, n), 2*n+1) for n in range(N+1)]
assert g == [1] + mul(mul(g, g), g)[:N]
q = [-c for c in mul(g, g)]
q[0] += 2
assert q[0] == 1
b = [1] + [0] * N
for n in range(1, N+1):
    b[n] = -sum(q[k] * b[n-k] for k in range(1, n+1))
assert mul(q, b) == [1] + [0] * N
a = mul(b, b)
assert mul(mul(q, q), a) == [1] + [0] * N
assert a[:len(SOURCE)] == SOURCE
print(f'Fidelity: all {len(SOURCE)} source coefficients (n=0..22) match.')
print('Initial coefficients:', ','.join(map(str, a[:9])))

# Independent stated source formula: Pell-weighted binomial sum, n>0.
pell = [0, 1]
for k in range(N+2):
    pell.append(2*pell[-1] + pell[-2])
for n in range(1, N+1):
    pell_sum = Fraction(sum(k * ((k+2)*pell[k+1] + (k+1)*pell[k+2])
                            * comb(3*n, n-k) for k in range(1, n+1)), 4*n)
    assert pell_sum == a[n], (n, pell_sum, a[n])
print(f'Fidelity: source Pell/binomial formula matches for n=1..{N}.')

def weights(n):
    return [10*n*(2*n-5), -511*n*n+1968*n-1337,
            2*(2461*n*n-12759*n+14814),
            4*(-5447*n*n+35358*n-55483),
            16*(2651*n*n-20323*n+38832),
            -2976*(3*n-11)*(3*n-13)]

boundary = [weights(5)[k]*a[5-k] for k in range(6)]
assert boundary == [2090000, -5164848, 4415488, -1546272, 223488, -23808]
assert sum(boundary) == -5952
print('Boundary n=5 summands:', ','.join(map(str, boundary)))
print('Boundary n=5 residual:', sum(boundary))
for n in range(6, N+1):
    assert sum(weights(n)[k]*a[n-k] for k in range(6)) == 0, n
    assert weights(n)[0] != 0
print(f'Recurrence: zero residual at every n=6..{N}.')

# Exact algebraic certificate; x and t are independent indeterminates.
x, t, n = s.symbols('x t n')
v, qt = 3-2*t, 2-t*t
A = 1/qt**2
tprime = t**4/v
Aprime = s.diff(A, t)*tprime
Asecond = s.diff(Aprime, t)*tprime
assert s.cancel(Aprime-4*t**5/(v*qt**3)) == 0
assert s.cancel(Asecond-4*t**8*(30-16*t+3*t*t-4*t**3)/(v**3*qt**4)) == 0
U0 = -8*x*(2976*x**4+88*x**3-784*x*x+215*x-15)
U1 = -2*(26784*x**5-7080*x**4-5352*x**3+2915*x*x-473*x+25)
U2 = -(27*x-4)*(8*x*x-8*x+1)*(124*x*x-54*x+5)
forcing = -20*x*x+208*x**3+320*x**4-5952*x**5
operator = U2*x*x*Asecond + x*(U1+U2)*Aprime + U0*A
assert s.cancel((operator-forcing).subs(x, (t-1)/t**3)) == 0
numerator, denominator = s.fraction(s.cancel(operator-forcing))
relation = x*t**3-t+1
quotient, remainder = s.div(s.Poly(numerator, x, t, domain=s.QQ),
                           s.Poly(relation, x, t, domain=s.QQ))
assert remainder.is_zero
assert s.expand(numerator - relation*quotient.as_expr()) == 0
assert denominator.subs({x: 0, t: 1}) != 0
print('Certificate: derivative formulas and rational substitution identity verified.')
print('Certificate: cleared numerator divisible by x*t^3-t+1, remainder=0.')
print('Certificate: quotient bidegree (x,t)=', quotient.degree(x), quotient.degree(t))

bracket = 5*t**9-52*t**7+52*t**6-80*t**5+160*t**4+1408*t**3-4464*t*t+4464*t-1488
structured = 5*t**9-52*(t-1)*t**6-80*(t-1)**2*t**3+1488*(t-1)**3
assert s.expand(bracket-structured) == 0
assert s.cancel(operator.subs(x, (t-1)/t**3) + 4*(t-1)**2/t**15*bracket) == 0
rows = [(s.expand(U2).coeff(x, k), s.expand(U1).coeff(x, k),
         s.expand(U0).coeff(x, k)) for k in range(6)]
assert rows == [(20, -50, 0), (-511, 946, 120), (4922, -5830, -1720),
                (-21788, 10704, 6272), (42416, 14160, -704),
                (-26784, -53568, -23808)]
for k, (u2, u1, u0) in enumerate(rows):
    assert s.expand(u2*(n-k)**2 + u1*(n-k) + u0 - weights(n)[k]) == 0
assert s.degree(forcing, x) == 5
assert s.diff(forcing, x, 6) == 0
# A nonzero leading coefficient remains after six differentiations.
assert not s.Poly(U2*x*x, x).is_zero
print('Euler coefficient rows:', rows)
print('Certificate: all six coefficient polynomials equal the conjectured weights.')
print('Certificate: forcing degree=5; six derivatives annihilate the forcing.')
print('Scope: finite coefficient checks and exact polynomial certificate; no Lean proof.')
print('ALL_CHECKS_OK')
