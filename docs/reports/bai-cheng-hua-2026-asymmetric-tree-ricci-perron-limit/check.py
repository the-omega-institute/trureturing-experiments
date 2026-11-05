#!/usr/bin/env python3
"""Exact checks for arXiv:2605.30949v1, Remark 5(2) and Section 5.2.

Only integers, Fraction and exact SymPy polynomial algebra are used. Spectral
intervals are certified by Sylvester's criterion on a symmetric Jacobi matrix
similar to the leaf-orbit quotient. Positive off-diagonal entries put the
Perron eigenvector in this quotient; omitted leaf-difference eigenvalues are
negative. Finite checks do not establish the universal convergence conjecture.
"""
from fractions import Fraction as F
from math import lcm
from pathlib import Path
import sympy as S

X = S.Symbol('X')

def ricci(p, q, length):
    edges = [(j, j + 1) for j in range(length)]
    edges += [(0, length + 1 + j) for j in range(p)]
    edges += [(length, length + 1 + p + j) for j in range(q)]
    deg = {}
    for e in edges:
        for v in e:
            deg[v] = deg.get(v, 0) + 1
    def entry(i, j):
        e, f = edges[i], edges[j]
        if i == j:
            return -sum((S.Rational(1, deg[v]) for v in e), S.Integer(0))
        return sum((S.Rational(1, deg[v]) for v in set(e) & set(f)), S.Integer(0))
    return S.Matrix(len(edges), len(edges), entry)

def quotient(p, q, length):
    d, e = p + 1, q + 1
    a = [-F(3, d)]
    if length == 1:
        a += [-F(1, d) - F(1, e)]
    else:
        a += [-F(1, d) - F(1, 2)] + [-F(1)] * (length - 2)
        a += [-F(1, e) - F(1, 2)]
    a += [-F(3, e)]
    b = [F(p, d*d)] + [F(1, 4)] * (length - 1) + [F(q, e*e)]
    scale = lcm(2, d, e)
    ai = [int(t*scale) for t in a]
    bi = [int(t*scale*scale) for t in b]
    assert all(F(i, scale) == t for i, t in zip(ai, a))
    assert all(F(i, scale*scale) == t for i, t in zip(bi, b))
    return scale, ai, bi

def above_top(x, data):
    scale, a, b = data
    num, den = x.numerator * scale, x.denominator
    older, prev = 1, num - den*a[0]
    if prev <= 0:
        return False
    for i in range(1, len(a)):
        current = (num - den*a[i])*prev - den*den*b[i-1]*older
        if current <= 0:
            return False
        older, prev = prev, current
    return True

def top_interval(p, q, length, bits=144):
    data = quotient(p, q, length)
    lo, hi = F(-2), F(2)
    assert not above_top(lo, data) and above_top(hi, data)
    for _ in range(bits):
        mid = (lo + hi)/2
        if above_top(mid, data):
            hi = mid
        else:
            lo = mid
    return lo, hi

def beta(s):
    return 2*X + 1 + S.Rational(2, s+1) - S.Rational(2*s, s+1)/((s+1)*X + 3)

def match_polynomial(s):
    # q + q^(-1) = 2(1+X). Retain the decaying root, beta in (0,1).
    b = beta(s)
    pol = S.Poly(S.cancel(b*b - 2*(1+X)*b + 1).as_numer_denom()[0], X)
    primitive = pol.primitive()[1]
    return primitive if primitive.LC() > 0 else -primitive

def limit_interval(s):
    pol = match_polynomial(s)
    roots = [(F(a), F(b)) for (a,b), multiplicity in pol.intervals(eps=S.Rational(1,2**155)) if a > 0]
    assert len(roots) == 1 and s > 3
    a,b = roots[0]
    assert 0 < beta(s).subs(X, S.Rational(a.numerator,a.denominator)) < 1
    assert 0 < beta(s).subs(X, S.Rational(b.numerator,b.denominator)) < 1
    return a,b

def fixed_decimal(x, digits=12):
    # Decimal display by integer truncation, not floating point.
    signed = x.numerator * 10**digits // x.denominator
    sign = '-' if signed < 0 else ''
    z = abs(signed)
    return f'{sign}{z//10**digits}.{z%10**digits:0{digits}d}'

def main():
    # Fidelity: Section 5.2's flat S_(2,5) and interface data.
    r = ricci(2, 5, 1)
    w = S.Matrix([3] + [1]*7)
    assert r*w == S.zeros(8,1)
    assert S.Rational(5-3, 6*3) == S.Rational(1,9)
    assert S.Rational(2-3, 3*3) == -S.Rational(1,9)
    assert beta(2).subs(X,0) == S.Rational(11,9)
    assert beta(5).subs(X,0) == S.Rational(7,9)
    assert beta(3).subs(X,0) == 1
    print('Fidelity: S_(2,5) length 1 has positive kernel vector (3,1,1,1,1,1,1,1); top = 0.')
    print('Fidelity: u=1/9, v=-1/9, Phi_A(0)=11/9, Phi_B(0)=7/9.')
    assert match_polynomial(5) == S.Poly(216*X**3 + 288*X**2 + 111*X - 2, X)
    pol = match_polynomial(5)
    assert pol.eval(S.Rational(1,60)) < 0 < pol.eval(S.Rational(1,50))
    assert all(c > 0 for c in S.Poly(S.diff(pol.as_expr(),X),X).all_coeffs())
    print('Matching polynomial: 216*X^3 + 288*X^2 + 111*X - 2; unique positive root in (1/60,1/50).')
    # Independently compare quotient characteristic polynomial with full matrix.
    for p,q in [(2,5),(4,6),(1,4)]:
        for length in range(1,5):
            scale,ai,bi=quotient(p,q,length)
            ds=[S.Rational(t,scale) for t in ai]
            bs=[S.Rational(t,scale*scale) for t in bi]
            older,prev=S.Integer(1),X-ds[0]
            for i in range(1,len(ds)):
                older,prev=prev,S.expand((X-ds[i])*prev-bs[i-1]*older)
            expected=prev*(X+1+S.Rational(2,p+1))**(p-1)*(X+1+S.Rational(2,q+1))**(q-1)
            assert S.expand(ricci(p,q,length).charpoly(X).as_expr()-expected)==0
    print('Fidelity cross-check: full rational Ricci characteristic polynomial = quotient times leaf factors (12 cases).')
    rows=['p,q,length,lambda_lower,lambda_upper,limit_lower,limit_upper']
    for p,q in [(2,5),(4,6),(1,4)]:
        limlo,limhi=limit_interval(max(p,q))
        previous=None
        for length in range(1,201):
            lo,hi=top_interval(p,q,length)
            rows.append(f'{p},{q},{length},{lo},{hi},{limlo},{limhi}')
            if length>=2:
                if (p,q)==(4,6):
                    assert hi < previous[0] and lo > limhi
                else:
                    assert lo > previous[1] and hi < limlo
            previous=(lo,hi)
            if length in [1,2,3,4,5,10,20,50,100,200]:
                print(f'S_({p},{q}) length={length}: certified [{fixed_decimal(lo)}, {fixed_decimal(hi+F(1,10**12))}] (exact endpoints in intervals.csv)')
        assert max(abs(lo-limhi),abs(hi-limlo)) < F(1,10**10)
        print(f'S_({p},{q}) length 200: |lambda-L| < 1/10^10; L in [{fixed_decimal(limlo)}, {fixed_decimal(limhi+F(1,10**12))}].')
    # The S_(5,5) source limit equation is literally the same beta_5 equation.
    assert S.cancel(beta(5) - (2*X+S.Rational(4,3)-S.Rational(5,9)/(2*X+1))) == 0
    # Source sanity check: S_(3,3) has a positive zero eigenvector for all these lengths.
    for length in range(1,201):
        scale,ai,bi=quotient(3,3,length)
        vals=[1]+[3]*length+[1]
        # nonsymmetric quotient equations (without square roots).
        for i in range(len(ai)):
            total = F(ai[i], scale)*vals[i]
            if i > 0:
                total += (F(3,4) if i == 1 else F(1,4) if i == len(ai)-1 else F(1,2))*vals[i-1]
            if i+1 < len(ai):
                total += (F(1,4) if i == 0 else F(3,4) if i == len(ai)-2 else F(1,2))*vals[i+1]
            assert total == 0
        assert not above_top(F(0), (scale,ai,bi))
        assert above_top(F(1,10**40), (scale,ai,bi))
    print('Fidelity: S_(3,3) top = 0 at every length 1..200 (positive vector and exact top enclosure).')
    Path(__file__).with_name('intervals.csv').write_text('\n'.join(rows)+'\n')
    print('Scope: 600 exact spectral enclosures; no universal theorem inferred from this finite check.')
    print('ALL_CHECKS_OK')

if __name__=='__main__':
    main()
