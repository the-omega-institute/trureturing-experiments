"""Yukawa mismatch required by exact three-channel protection at the nominal 1 TeV input."""
from flint import arb, arb_mat, ctx


def dot(a, b):
    return (a.transpose() * b)[0, 0]


def minimum_skew(up, charm, top, x, v, w):
    def column(a, b):
        return [charm*a[0, 0], top*b[0, 0], top*b[1, 0] - charm*a[2, 0]]
    plus = arb_mat([column(v, w), column(w, -v)]).transpose()
    minus = arb_mat([column(v, -w), column(w, v)]).transpose()
    target = up * arb_mat([[x[1, 0]], [x[2, 0]], [0]])

    def vectors(t):
        return arb_mat([[t.cos()], [t.sin()]]), arb_mat([[-t.sin()], [t.cos()]])
    def slope(t):
        h, dh = vectors(t)
        return dot(minus*dh, minus*h - target)
    lo, hi = arb('.17456'), arb('.17458')
    assert slope(lo) < 0 < slope(hi)
    for _ in range(70):
        mid = (lo + hi) / 2
        value = slope(mid)
        if value > 0:
            hi = mid
        else:
            assert value < 0
            lo = mid
    theta = lo.union(hi)
    h, dh = vectors(theta)
    ell = dot(target, minus*h)
    assert ell > 0
    def reduced(m):
        v = m.transpose()*target
        return m.transpose()*m - v*v.transpose()/ell
    s = reduced(minus)
    lam = dot(h, s*h)
    assert dot(dh, s*dh) > lam
    other = reduced(plus)
    low_eigenvalue = (other.trace() - ((other[0, 0] - other[1, 1])**2 + 4*other[0, 1]**2).sqrt())/2
    assert low_eigenvalue > lam
    optimum = dot(minus*h - target, minus*h - target).sqrt()
    assert arb('.86080564118205') < optimum < arb('.86080564118207')
    return optimum, theta


def main():
    ctx.prec = 192
    up, charm, top = map(arb, ('.00000615', '.00311', '.8616'))
    down, strange, bottom = map(arb, ('.0000135', '.000268', '.01401'))
    electron, muon, tau = map(arb, ('.0000028227', '.00059465', '.010101'))
    t13, t23 = arb('.00377'), arb('.04275')
    p, q, r = t13.sin(), t23.sin()*t13.cos(), t23.cos()*t13.cos()
    t = t13.cos()
    lower = ((charm*q)**2 + (top*r)**2).sqrt() - up*t

    # Complete the light right-handed direction to a real orthonormal mass frame.
    x = arb_mat([[p], [q], [r]])
    v = arb_mat([[t], [-p*q/t], [-p*r/t]])
    w = arb_mat([[0], [-r/t], [q/t]])
    optimum, theta = minimum_skew(up, charm, top, x, v, w)
    assert arb('.86080054') < lower < optimum
    assert optimum / top > arb('.9990')

    protected_scale = (bottom + down) / (muon + tau)
    protected = bottom - protected_scale*muon
    assert abs(strange - protected_scale*electron) < protected
    assert arb('.01323033029') < protected < arb('.01323033030')
    assert protected/bottom > arb('.9443')
    free_scale = (bottom + strange) / (muon + tau)
    free = bottom - free_scale*tau
    assert abs(down - free_scale*electron) < free
    assert arb('.00052581923') < free < arb('.00052581924')
    print('Down-lepton mismatch, protected', protected.str(16), 'scale', protected_scale.str(16))
    print('Down-lepton mismatch, free', free.str(16), 'scale', free_scale.str(16))
    print('Up transpose mismatch lower', lower.str(16), 'sharp minimum', optimum.str(16), 'rotation', theta.str(16))


if __name__ == '__main__':
    main()
