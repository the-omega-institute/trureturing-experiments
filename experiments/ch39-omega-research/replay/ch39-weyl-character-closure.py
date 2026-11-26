"""Two-cell minimax spread for a central weight on the five native Weyl multiplets."""
from fractions import Fraction as Q
from flint import arb, ctx

ctx.prec = 192


def responses(powers):
    a, u, d, l, e = powers
    return (3*(a+4*u+d+l+e)/5, 9*a+l, 6*a+3*(u+d)/2)


def minimax():
    root = arb(71).sqrt()
    s0 = 6 + 2*root
    R = (742-22*root)/445
    t = -arb(17)/178 + 671*root/6319
    assert 0 < t < 1 and arb('1.25084203368') < R < arb('1.25084203369')

    # H = t A + (1-t) B - R C, with Q amplitude one, in variables r,s.
    L = arb(9)/4 * (1 - 2*t/5)
    rr, ss, rs = 39*t/5 + L - 3*R/2, 3*t/5 + L - 3*R/2, 2*L
    r1, s1 = -216*t/5 - 4*L, -4*L
    c = 87*t + 9*(1-t) + 4*L - 6*R
    # For r <= 2 minimize over all s. The minimum is a concave quadratic in r.
    assert ss > 0 and rr - rs**2/(4*ss) < 0
    assert 4*rr + 2*r1 + c > 0
    # For r >= 2, s >= 1, H increases with s. H(r,1) is strictly positive.
    assert rr > 0 and 4*rr*(ss+s1+c) - (rs+r1)**2 > 0
    # For 2 <= r <= 4, 0 <= s <= 1, A - R B is positive.
    assert 3*(1+arb(576)/13)/5 - 9*R - (R-arb(3)/5)*arb(81)/4 > 0
    vals = responses([arb(1), arb(0), s0**2, (3*s0/2-3)**2, arb(144)])
    assert (vals[0]-vals[1]).contains(0) and (vals[0]/vals[2]-R).contains(0)
    print('Sharp two-cell max/min ratio', R.str(16), 'boundary d/Q', s0.str(14))


def ordinary_trace():
    # The factors d_i I_ia carry a common conjugate-pair multiplier.
    dims = [6, 3, 3, 2, 1]
    for cells in (1, 2, 3):
        result = responses([Q(1, d)**cells for d in dims])
        if cells == 1:
            assert result[0] == result[1] == result[2] == 2
        elif cells == 2:
            assert result == (Q(11,10), Q(1,2), Q(1,2))
        else:
            assert result[0] > result[1] > result[2]
    print('Unweighted Weyl ledger, two-cell ratio 11:5:5.')


def extra_characters():
    # Scaled differences 48*(weak-color), 240*(hypercharge-color).
    def defect(p, q, n, y):
        d = (p+1)*(q+1)*(p+q+2)*n//2
        casimir = p*p + q*q + p*q + 3*(p+q)
        orbit = 1 if p == q and y == 0 else 2
        return (orbit*d*d*(2*(n*n-1)-casimir), orbit*d*d*(2*y*y-5*casimir))
    extra = [(0,0,4,-6), (0,1,2,-5), (1,1,1,-3)]
    assert tuple(map(sum, zip(*(defect(*rep) for rep in extra)))) == (-48,1008)
    print('Three additional character pairs restore all-cell equality at equal eigenvalues.')


if __name__ == '__main__':
    ordinary_trace()
    minimax()
    extra_characters()
