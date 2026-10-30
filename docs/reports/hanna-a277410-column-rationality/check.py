#!/usr/bin/env python3
"""Exact finite checks for Hanna's A277410 column conjecture (OEIS revision 64).

Two independent constructions: truncated source fixed-point iteration in Q[y][t],
and columnwise exponential-polynomial Taylor recurrence in Q[t, exp(t)].
All arithmetic carrying evidence is exact. These checks are not an all-k proof.
"""
from math import factorial
from pathlib import Path
import sys
import sympy as s

# Verbatim numerical examples from https://oeis.org/A277410, EXAMPLE.
SOURCE_ROWS = [[1], [1, 0], [1, 3, 0], [1, 13, 15, 0], [1, 38, 165, 105, 0], [1, 94, 1033, 2310, 945, 0], [1, 213, 4953, 26229, 36330, 10395, 0], [1, 459, 20370, 213511, 674520, 640710, 135135, 0], [1, 960, 76056, 1421225, 8559675, 18127935, 12588345, 2027025, 0], [1, 1972, 266334, 8283234, 85654979, 337805535, 515903850, 273544425, 34459425, 0], [1, 4007, 892542, 44013478, 729292193, 4822487682, 13506364410, 15631793100, 6529047525, 654729075, 0], [1, 8089, 2900353, 218797958, 5531376285, 57226590953, 264482764305, 555756298020, 505173143475, 170116046100, 13749310575, 0]]
SOURCE_P = ['1', '3*x - 2*x^2', '15*x - 45*x^2 - 2*x^3 + 106*x^4 - 92*x^5 + 24*x^6', '105*x - 840*x^2 + 504*x^3 + 16321*x^4 - 75880*x^5 + 154483*x^6 - 152077*x^7 + 39208*x^8 + 59000*x^9 - 60336*x^10 + 23328*x^11 - 3456*x^12', '945*x - 15645*x^2 + 32445*x^3 + 1255770*x^4 - 15120061*x^5 + 86803308*x^6 - 291640845*x^7 + 529758178*x^8 - 50236668*x^9 - 2553002523*x^10 + 7695202852*x^11 - 12713196156*x^12 + 13351222596*x^13 - 8752472980*x^14 + 2871967920*x^15 + 387984096*x^16 - 884504448*x^17 + 427064832*x^18 - 100694016*x^19 + 9953280*x^20', '10395*x - 305235*x^2 + 1299375*x^3 + 77300220*x^4 - 1834009998*x^5 + 21447595316*x^6 - 156933684108*x^7 + 721294719700*x^8 - 1490891586137*x^9 - 5868653004882*x^10 + 70213320019895*x^11 - 359261247450016*x^12 + 1234731543184308*x^13 - 3081038591203028*x^14 + 5553265322783926*x^15 - 6518085613542516*x^16 + 2256970375232288*x^17 + 9498116639867573*x^18 - 25485484994020128*x^19 + 37162639109810884*x^20 - 37419816866322296*x^21 + 27200926921683600*x^22 - 14055671260790656*x^23 + 4698364855901568*x^24 - 583485067952640*x^25 - 341605998065664*x^26 + 237336648708096*x^27 - 72380729917440*x^28 + 11910492979200*x^29 - 859963392000*x^30']
t, y, z = s.symbols('t y z')
N, K = 12, 6
ZP = s.Poly(0, t, domain=s.QQ)
YP = s.Poly(0, y, domain=s.QQ)

def poly(value):
    return s.Poly(value, t, domain=s.QQ)

def clean(a):
    return {j: p for j, p in a.items() if not p.is_zero}

def add(a, b):
    out = dict(a)
    for j, p in b.items():
        out[j] = out.get(j, ZP) + p
    return clean(out)

def scale(a, c):
    return clean({j: p.mul_ground(c) for j, p in a.items()})

def mul(a, b):
    out = {}
    for j, p in a.items():
        for l, q in b.items():
            out[j+l] = out.get(j+l, ZP) + p*q
    return clean(out)

def deriv(a, count=1):
    for _ in range(count):
        a = clean({j: p.diff() + p.mul_ground(j) for j, p in a.items()})
    return a

def integrate(a):
    out = {}
    for j, p in a.items():
        if j == 0:
            out[j] = p.integrate()
        else:
            q, d = ZP, p
            for i in range(int(p.degree())+1):
                q += d.mul_ground(s.Rational((-1)**i, j**(i+1)))
                d = d.diff()
            out[j] = q
    initial = sum(p.eval(0) for p in out.values())
    out[0] = out.get(0, ZP) - poly(initial)
    return clean(out)

def compositions(n, size):
    if size == 1:
        yield (n,)
    else:
        for a in range(n+1):
            for rest in compositions(n-a, size-1):
                yield (a,) + rest

def ode_solve(forcing):
    out = {}
    for j, p in forcing.items():
        c = j-1
        if c == 0:
            out[j] = p.integrate()
        else:
            q, d = ZP, p
            for i in range(int(p.degree())+1):
                q += d.mul_ground(s.Rational((-1)**i, c**(i+1)))
                d = d.diff()
            out[j] = q
    out[1] = out.get(1, ZP) - poly(sum(p.eval(0) for p in out.values()))
    return clean(out)

def factorial_coeff(a, n):
    value = 0
    for j, p in a.items():
        for (d,), c in p.terms():
            if d <= n:
                value += c * (factorial(n)//factorial(n-d)) * j**(n-d)
    return s.Rational(value)

def tmul(a, b):
    out = [YP for _ in range(N+1)]
    for i, p in enumerate(a):
        if not p.is_zero:
            for j in range(N+1-i):
                if not b[j].is_zero:
                    out[i+j] += p*b[j]
    return out

def tcompose(a, b):
    out = [YP for _ in range(N+1)]
    for p in reversed(a):
        out = tmul(out, b)
        out[0] += p
    return out

def tintegral(a):
    return [YP] + [a[n].mul_ground(s.Rational(1,n+1)) for n in range(N)]

def check_source():
    x = [YP, s.Poly(1,y,domain=s.QQ)] + [YP]*(N-1)
    a = list(x)
    for _ in range(N+1):
        inner = [p.mul_ground(1) for p in x]
        for n in range(N+1):
            inner[n] += (a[n]-x[n])*s.Poly(y,y,domain=s.QQ)
        b = tcompose(tintegral(a),inner)
        b[1] += s.Poly(1,y,domain=s.QQ)
        a = b
    # Check the original NAME equation, rather than only its fixed-point form.
    g = tintegral(a)
    left_inner = [x[n]-g[n]*s.Poly(y,y,domain=s.QQ) for n in range(N+1)]
    left = tcompose(a,left_inner)
    right = [x[n]+g[n]*s.Poly(1-y,y,domain=s.QQ) for n in range(N+1)]
    assert left == right
    count = 0
    for n,row in enumerate(SOURCE_ROWS,1):
        for k,expected in enumerate(row):
            assert a[n].nth(k)*factorial(n) == expected, (n,k)
            count += 1
        assert a[n].degree() < n
    print(f'FIDELITY: {count} source triangle entries (rows 1..{N}) match')
    print(f'SOURCE_EQUATIONS: original and fixed-point equations hold modulo t^{N+1}')
    return a

def check_columns(a):
    fs = [{0:poly(-1),1:poly(1)}]
    hs = [add(fs[0],{0:poly(-t)})]
    for k in range(1,K+1):
        remainder = {}
        for r in range(1,k+1):
            for indices in compositions(k-r,r+1):
                term = deriv(fs[indices[0]],r-1)
                for ai in indices[1:]:
                    term = mul(term,hs[ai])
                remainder = add(remainder,scale(term,s.Rational(1,factorial(r))))
        assert remainder.get(0,ZP) == poly(t+k)
        for j,p in remainder.items():
            if j > 0:
                assert j <= k+1 and p.degree() <= 2*(k+1-j)
            if j == 1:
                assert p.degree() <= 2*k-1
        f = ode_solve(deriv(remainder))
        assert add(deriv(f),scale(f,-1)) == deriv(remainder)
        assert add(f,scale(integrate(f),-1)) == remainder
        assert sum(p.eval(0) for p in f.values()) == 0
        fs.append(f)
        hs.append(f)
    for k,f in enumerate(fs):
        assert f.get(0,ZP) == poly(-1)
        assert all(j <= k+1 and p.degree() <= 2*(k+1-j)
                   for j,p in f.items() if j > 0)
        for n in range(N+1):
            assert factorial_coeff(f,n) == a[n].nth(k)*factorial(n), (k,n)
        D = s.Poly(s.prod((1-j*z)**(2*(k+1-j)+1) for j in range(1,k+2)),z,domain=s.QQ)
        numerator = -D
        for j,p in f.items():
            if j == 0:
                continue
            for (d,),c in p.terms():
                divisor = s.Poly((1-j*z)**(d+1),z,domain=s.QQ)
                numerator += D.exquo(divisor)*s.Poly(c*factorial(d)*z**d,z,domain=s.QQ)
        assert numerator.degree() == (k+1)**2
        assert numerator.LC() == -D.LC()
        P = numerator.exquo(s.Poly(z**(k+1),z,domain=s.QQ))
        assert P.degree() == k*(k+1)
        if k < len(SOURCE_P):
            expected = s.Poly(s.sympify(SOURCE_P[k].replace('^','**'),locals={'x':z}),z,domain=s.QQ)
            assert P == expected, k
        if k > 0:
            assert all(factorial_coeff(f,n)==0 for n in range(k+2))
        # Independently multiply D by the factorial coefficients, through 20
        # coefficients beyond the claimed numerator degree.
        for n in range((k+1)**2+21):
            value = sum(D.nth(i)*factorial_coeff(f,n-i) for i in range(min(n,int(D.degree()))+1))
            assert value == numerator.nth(n), (k,n)
        print(f'COLUMN k={k}: rates/degrees, ODE, initial zeros, rational identity OK; '
              f'deg(D*C)={numerator.degree()}, deg(P)={P.degree()}, LC={numerator.LC()}')
    print('SOURCE_NUMERATORS: P_0,...,P_5 reproduced exactly')
    print('SCOPE: finite exact checks for k=0,...,6; no uniform proof claimed')

if __name__ == '__main__':
    print(f'Python {sys.version.split()[0]}; SymPy {s.__version__}; exact rational arithmetic')
    check_columns(check_source())
    print('ALL_CHECKS_OK')
