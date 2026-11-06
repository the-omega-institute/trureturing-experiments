#!/usr/bin/env python3
"""Exact checks of Puddu, arXiv:2611.06638v1, Conjecture 4.1.

Run: python3 check.py
Requires SymPy. All polynomial and root calculations are over QQ, with rational
root-isolating intervals. Local quotient dimensions are recomputed by Groebner
bases, without replacing the equations by the proposed monomial ideal.
"""
from collections import defaultdict
from hashlib import sha256
from pathlib import Path
import sympy as s
import subprocess
import sys

sys.stdout.reconfigure(line_buffering=True)

x, w, t, a, b = s.symbols('x w t a b')


def derivative(p, n):
    q = defaultdict(int)
    for (i, j), c in p.items():
        if i:
            q[i - 1, j + n] += c * i
        if j:
            q[i + n, j - 1] += c * j
    return dict(q)


def expression(p):
    return sum(c * a**i * b**j for (i, j), c in p.items())


def primitive(p, n, j):
    u, v = (-j) % (n + 1), (1 - j) % (n + 1)
    d = 1 + (n - 1) * j
    h, rem = divmod(d - u - v, n + 1)
    assert rem == 0 and h >= 0
    assert min(q for i, q in p) == u
    assert min(i for i, q in p) == v
    assert all(i + q == d and q % (n + 1) == u and i % (n + 1) == v
               and c > 0 for (i, q), c in p.items())
    P = s.Poly.from_dict({((q - u) // (n + 1),): c for (i, q), c in p.items()}, t, domain=s.QQ)
    assert P.degree() == h and P.nth(0) > 0 and P.LC() > 0
    Q = s.Poly(sum(P.nth(i) * t**(h - i) for i in range(h + 1)), t, domain=s.QQ)
    return u, v, h, P, Q


def alternating(P, Q):
    if P.degree() == 0:
        return
    assert s.gcd(P, Q).degree() == 0
    eps = s.Rational(1, 10**12)
    while True:
        intervals = []
        for label, p in enumerate((P, Q)):
            roots = p.intervals(eps=eps)
            assert len(roots) == p.degree()
            for (lo, hi), mult in roots:
                assert mult == 1 and hi < 0
                intervals.append((lo, hi, label))
        intervals.sort()
        if all(intervals[i][1] < intervals[i + 1][0] for i in range(len(intervals) - 1)):
            break
        eps /= 10**6
    assert all(intervals[i][2] != intervals[i + 1][2] for i in range(len(intervals) - 1))


def chart_equation(jets, n, m, which, swap, limit):
    d = 1 + (n - 1) * (m + 1)
    # A rational specialization of the normalized hyperplanes A_top+U, B_top+V.
    terms = defaultdict(int)
    top = jets[m + 1]
    for (i, j), c in top.items():
        i, j = (j, i) if swap else (i, j)
        terms[j, 0] += c
    for ell in range(m + 1):
        for coord in range(2):
            p = jets_all[coord][ell]
            coefficient = 2 + (which + 1) * (2 * ell + coord + 1)
            for (i, j), c in p.items():
                i, j = (j, i) if swap else (i, j)
                z = d - (i + j)
                if j < limit and z < limit:
                    terms[j, z] += coefficient * c
    if d < limit:
        terms[0, d] += 7 + which
    return s.Poly.from_dict({(i, j): c for (i, j), c in terms.items()
                            if i < limit and j < limit}, x, w, domain=s.QQ)


def local_length(n, m, swap, k, H):
    lengths = []
    for L in (k + H + 1, k + H + 2):
        f = chart_equation(jets_all[0], n, m, 0, swap, L)
        g = chart_equation(jets_all[1], n, m, 1, swap, L)
        gb = s.groebner([f.as_expr(), g.as_expr(), x**L, w**L], x, w, order='grevlex', domain=s.QQ)
        leading = [p.LM(order=gb.order).exponents for p in gb.polys]
        count = sum(not any(i >= u and j >= v for u, v in leading)
                    for i in range(L) for j in range(L))
        # This tests the actual quotient, not just its dimension.
        assert {s.expand(p.as_expr()) for p in gb.polys} == {x**k, w**H}, (n, m, swap, gb)
        lengths.append(count)
    assert lengths == [k * H, k * H]
    return lengths[0]


def global_length(n, m):
    """Independent global QQ quotient dimension for an explicit affine section.

    The section is not required to be reduced. Its two local unit conditions
    at infinity hold whenever r is interior. This is an intersection-length
    check; the generic distinct-point interpretation still uses transversality.
    """
    r = m % (n + 1)
    f, g = expression(jets_all[0][m + 1]), expression(jets_all[1][m + 1])
    if 0 < r < n:
        f += expression(jets_all[0][m + 1 - r])
        g += 2 * expression(jets_all[1][m + 1 - r])
    source = ('ring R=0,(a,b),dp;\npoly f=' + str(f).replace('**', '^')
              + ';\npoly g=' + str(g).replace('**', '^')
              + ';\nideal I=std(ideal(f,g));\nprint(vdim(I));\nquit;\n')
    process = subprocess.run(['Singular', '-q'], input=source, text=True,
                             capture_output=True, check=True)
    assert not process.stderr, process.stderr
    answer = process.stdout.strip()
    assert answer.isdecimal(), (n, m, answer)
    return int(answer)


print('Python exact algebra; SymPy', s.__version__)
print('check.py sha256', sha256(Path(__file__).read_bytes()).hexdigest())
print('Global Groebner backend: Singular, ring characteristic 0, coefficient field QQ')
print('Fidelity: source Section 4 first three jet pairs and cubic m=1,2 gcds')
for n in range(2, 6):
    A, B = [{(1, 0): 1}], [{(0, 1): 1}]
    for j in range(13):
        A.append(derivative(A[-1], n))
        B.append(derivative(B[-1], n))
    jets_all = A, B
    assert expression(A[1]) == b**n
    assert s.expand(expression(A[2]) - n * a**n * b**(n - 1)) == 0
    assert s.expand(expression(A[3]) - n * a**(n - 1) * b**(n - 2)
                    * ((n - 1) * a**(n + 1) + n * b**(n + 1))) == 0
    for j in range(14):
        u, v, h, P, Q = primitive(A[j], n, j)
        b_dehom = s.Poly(sum(c * x**q for (i, q), c in B[j].items()), x)
        assert b_dehom == s.Poly(x**v * Q.as_expr().subs(t, x**(n + 1)), x)
        alternating(P, Q)
    for m in range(13):
        r = m % (n + 1)
        k = n - r if 0 < r < n else 0
        gcd = s.Poly(s.gcd(expression(A[m + 1]), expression(B[m + 1])), a, b, domain=s.QQ).monic()
        assert gcd == s.Poly((a * b)**k, a, b, domain=s.QQ)
        la = lb = 0
        if k:
            H = (n - 1) * r
            la = local_length(n, m, False, k, H)
            lb = local_length(n, m, True, k, H)
        d = (n - 1) * m + n
        degree = d*d - la - lb
        expected = d*d - 2*(n - 1)*r*(n - r)
        assert degree == expected
        global_dim = global_length(n, m)
        assert global_dim == expected, (n, m, global_dim, expected)
        if n == 2:
            assert degree == (m + 2)**2 - (2 if m % 3 == 1 else 0)
        if n == 3:
            assert degree == (2*m + 3)**2 - (8 if m % 4 in (1, 2) else 0)
        print(f'n={n} m={m} gcd=(ab)^{k} local_lengths={la},{lb} Bezout_degree={degree} global_QQ_length={global_dim}')
print('52 gcd cases; rational root isolation through j=13; two rational local Groebner quotients per point')
print('52 independent global QQ Groebner quotient dimensions agree. Specialized sections can be nonreduced.')
print('Generic distinct-point degree follows from the stated Bezout/transversality bridge; it is not a direct generic-point enumeration.')
print('ALL_CHECKS_OK')
