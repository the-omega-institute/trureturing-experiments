"""Exact triangle averages for one normalized Weyl generation, using rational arithmetic."""
from fractions import Fraction as F

DIMS = (6, 3, 3, 2, 1)
EPS = F(1, 100)
X = tuple(F(i, 1380) for i in (11, 41, 835, 1281, 9))
Z = tuple(F(i, 111) for i in (8, 2, 40, 39, 90))


def product(x, y):
    return tuple(i*j for i, j in zip(x, y, strict=True))


def response(x):
    q, u, d, l, e = x
    return 3*(q+4*u+d+l+e)/5, 9*q+l, 6*q+3*(u+d)/2


def mass(x):
    return sum(2*d*d*i for d, i in zip(DIMS, x, strict=True))


def denominator(x, y):
    return (1-EPS*mass(x))*(1-EPS*mass(y))+EPS**2*mass(product(x, y))


def coefficients(x, y):
    # Raw coefficients on the trivial mode and five real normalized characters.
    pair = product(x, y)
    return ((1-EPS*mass(x))*(1-EPS*mass(y)),
            *(2*EPS**2*d*d*v for d, v in zip(DIMS, pair, strict=True)))


def edge(x, y):
    scale = 2*EPS**2/denominator(x, y)
    return tuple(scale*k for k in response(product(x, y)))


def third_vertex(h):
    return tuple(z+h*(z-x) for x, z in zip(X, Z, strict=True))


def averages(y):
    xx, xy = edge(X, X), edge(X, y)
    mean = tuple((a+2*b)/3 for a, b in zip(xx, xy, strict=True))
    dxx, dxy = denominator(X, X), denominator(X, y)
    from_mean_weight = tuple((dxx*a+2*dxy*b)/(dxx+2*dxy) for a, b in zip(xx, xy, strict=True))
    return mean, from_mean_weight


def verify():
    dxx, dxz = denominator(X, X), denominator(X, Z)
    h = dxz/(3*dxx-dxz)
    y = third_vertex(h)
    other = third_vertex(F(1, 2))
    for table in (X, Z, y, other):
        assert min(table) > 0 and EPS*mass(table) < F(1, 2)
        assert response(table) == (1,)*3
    assert len(set(edge(X, Z))) == 1
    assert h == F(8793328720, 16665497969)
    assert h == denominator(X, y)/(2*dxx)
    mean, different = averages(y)
    assert mean == (F(5639, 57347796),)*3 and len(set(different)) > 1
    xx, xy = edge(X, X), edge(X, y)
    assert F(1565, 1000) < max(xx)/min(xx) < F(1566, 1000)
    assert F(4209, 1000) < max(xy)/min(xy) < F(4210, 1000)
    different, from_mean_weight = averages(other)
    assert len(set(from_mean_weight)) == 1 and len(set(different)) > 1
    for i, j, k, n in zip(coefficients(X, X), coefficients(X, y), coefficients(X, Z),
                          coefficients(X, other), strict=True):
        assert (i/dxx+2*j/denominator(X, y))/3 == k/dxz
        assert (i+2*n)/3 == k
    fourth = tuple(2*(a-b)**2/3 for a, b in zip(xx, xy, strict=True))
    for value, lo, hi in zip(fourth, (249, 417, 75), (250, 418, 76), strict=True):
        assert lo*F(1, 10**10) < value < hi*F(1, 10**10)
    print('Common single-cell stiffness', 2*EPS)
    print('Mean of three edge stiffnesses', mean[0], 'spreads', max(xx)/min(xx), max(xy)/min(xy))
    print('Stiffness of the mean boundary weight with its own table', from_mean_weight[0])


if __name__ == '__main__':
    verify()
