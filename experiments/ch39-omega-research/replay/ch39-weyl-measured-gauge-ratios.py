"""Exact elimination for the first native Weyl blocking length matching the gauge ratios."""
from pathlib import Path
from runpy import run_path
from sympy import QQ, Poly, Rational as Q, groebner, symbols
from flint import arb, ctx

ctx.prec = 768
responses = run_path(str(Path(__file__).with_name('ch39-weyl-character-closure.py')))['responses']
ctx.prec = 768
s, r = symbols('s r')
g1, g2, g3 = map(Q, ('.467453', '.63811', '1.0583'))
target = [(g3/g1)**2, (g3/g2)**2, Q(1)]
beta = [Q(41,10), -Q(19,6), -Q(7)]


def ball(x):
    return arb(int(x.p))/int(x.q)


def branches(cells):
    ell, e = 3*(r+s)/2-3, 12-3*r
    a, b, c = responses([Q(1), r**cells, s**cells, ell**cells, e**cells])
    basis = groebner([a - target[0]*c, b - target[1]*c], s, r, domain=QQ)
    linear, univariate = basis.polys
    first = Poly(linear.as_expr(), s)
    assert first.degree() == 1 and first.LC() == 1 and univariate.degree(s) == 0
    profile = Poly(-first.TC(), r, domain=QQ)
    roots = Poly(univariate.as_expr(), r, domain=QQ)
    assert roots.degree() == cells**2 and roots.gcd(roots.diff()).degree() == 0
    found = []
    for (lo, hi), multiplicity in roots.intervals(eps=Q(1,10**75), inf=0, sup=4):
        assert multiplicity == 1
        rr = ball(lo).union(ball(hi))
        ss = arb(0)
        for coefficient in profile.all_coeffs():
            ss = rr*ss + ball(coefficient)
        if ss < 0:
            continue
        assert rr > 0 and rr < 4 and ss > 0 and 3*(rr+ss)/2 > 3
        found.append((rr, ss))
    return found


def running_plane(rr, ss, end):
    normal = [target[(i+1)%3]*beta[(i+2)%3] - target[(i+2)%3]*beta[(i+1)%3] for i in range(3)]
    normal = [ball(v/normal[0]) for v in normal]
    x = [arb(1), rr, ss, 3*(rr+ss)/2-3, 12-3*rr]
    assert min(rr.lower(), ss.lower()) > 1 and x[-1] > x[-2] > rr and x[-2] > ss
    weights = []
    for j in range(5):
        row = responses([Q(int(i == j)) for i in range(5)])
        weights.append(sum(t*ball(z) for t, z in zip(normal, row, strict=True)))
    assert [1 if v > 0 else -1 if v < 0 else 0 for v in weights] == [-1,1,1,-1,1]
    def value(n):
        return sum(a*b**n for a,b in zip(weights,x,strict=True))
    assert value(0) < 0 < value(1)
    assert value(3) > 0 > value(5)
    assert value(end) < 0 < value(end+1)


def main():
    assert not branches(2) and not branches(3)
    roots = branches(4)
    assert len(roots) == 2
    for (rr, ss), end in zip(roots, (6, 10), strict=True):
        running_plane(rr, ss, end)
        print('Four-cell solution u/Q, d/Q:', rr.str(14), ss.str(14),
              'last crossing between', end, 'and', end+1)


if __name__ == '__main__':
    main()
