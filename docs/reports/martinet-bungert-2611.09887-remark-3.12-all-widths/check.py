#!/usr/bin/env python3
"""Exact rational checks for Martinet--Bungert Remark 3.12.

Polygon integration certifies finite test instances, not the universal theorem.
Only Python's standard library is needed; no floating point is used.
"""
from fractions import Fraction as F
from itertools import combinations, product

ZERO = (F(0), F(0), F(0))
Y = (F(0), F(1), F(0))


def add(a, b):
    return tuple(x + y for x, y in zip(a, b))


def scale(c, a):
    return tuple(c * x for x in a)


def value(a, p):
    return a[0] * p[0] + a[1] * p[1] + a[2]


def signed_area(poly):
    return sum((a[0] * b[1] - b[0] * a[1]
                for a, b in zip(poly, poly[1:] + poly[:1])), F(0)) / 2


def clip(poly, a):
    """Intersect a convex CCW polygon with the closed half-plane a >= 0."""
    out = []
    for p, q in zip(poly, poly[1:] + poly[:1]):
        vp, vq = value(a, p), value(a, q)
        if vp >= 0:
            out.append(p)
        if (vp < 0 < vq) or (vq < 0 < vp):
            t = vp / (vp - vq)
            out.append(tuple(p[i] + t * (q[i] - p[i]) for i in range(2)))
    out = [p for i, p in enumerate(out) if not i or p != out[i - 1]]
    if len(out) > 1 and out[0] == out[-1]:
        out.pop()
    return out if len(out) >= 3 and signed_area(out) > 0 else []


def split(poly, a):
    vals = [value(a, p) for p in poly]
    if min(vals) >= 0:
        return [(poly, a)]
    if max(vals) <= 0:
        return [(poly, ZERO)]
    return [(q, b) for q, b in [(clip(poly, a), a),
                               (clip(poly, scale(-1, a)), ZERO)] if q]


def constants(n, L, source=False):
    if source:
        assert n == 4
        return F(1), F(1, 100), F(1), F(16)
    eps = F(1, 4 * n * n)
    b = eps / (n + 1)
    delta = eps / (16 * n * (n + 1))
    return b, delta, eps, F(4 * (L + 1))


def hidden(n, eta):
    return [(F(1), 2 * eta, F(-1)), (F(1), eta, F(-2))] + [
        (F(-1), (eta if i == n - 1 else 2 * eta if i == n else F(0)), F(i))
        for i in range(3, n + 1)]


def weights(n, L, source=False):
    b, delta, eps, lam = constants(n, L, source)
    left = [b, lam] + [F(0)] * (n - 2)
    right = [F(0), F(0), lam] + [2 * eps] * (n - 4) + [eps]
    return left, right, delta


def affine_sum(cs, hs):
    out = ZERO
    for c, h in zip(cs, hs):
        out = add(out, scale(c, h))
    return out


def outer(n, L, hs, eta, source=False):
    wl, wr, delta = weights(n, L, source)
    kl, kr = affine_sum(wl, hs), affine_sum(wr, hs)
    if eta:
        zr = scale(1 / eta, add(add(hs[0], scale(-1, hs[1])), (0, 0, -1)))
        zl = scale(1 / eta, add(add(hs[-1], scale(-1, hs[-2])), (0, 0, -1)))
    else:
        zl = zr = Y
    return [add(add(z, scale(-1, k)), (0, 0, -j * delta))
            for z, k in [(zl, kl), (zr, kr)] for j in range(1, 4)]


def direct(n, L, eta, p, source=False):
    hs = []
    for h in hidden(n, eta):
        hs.append(h if value(h, p) > 0 else ZERO)
    active = [a if value(a, p) > 0 else ZERO for a in outer(n, L, hs, eta, source)]
    out = affine_sum([1] * 6, active)
    return value(out, p), out[:2]


def exact_partition(n, L, eta, source=False):
    rectangle = [(F(-L), F(-L)), (F(n + 1 + L), F(-L)),
                 (F(n + 1 + L), F(L)), (F(-L), F(L))]
    cells = [(rectangle, [])]
    for h in hidden(n, F(0)) + hidden(n, eta):
        cells = [(q, hs + [ah]) for poly, hs in cells for q, ah in split(poly, h)]
    out = []
    for poly, hs in cells:
        states = [(poly, ZERO, ZERO)]
        args = outer(n, L, hs[:n], F(0), source) + outer(n, L, hs[n:], eta, source)
        for i, a in enumerate(args):
            states = [(q, add(f, ah) if i < 6 else f,
                       u if i < 6 else add(u, ah))
                      for p, f, u in states for q, ah in split(p, a)]
        out.extend(states)
    assert sum((signed_area(p) for p, _, _ in out), F(0)) == 2 * L * (n + 1 + 2 * L)
    return out


def integral_square(poly, a):
    out = F(0)
    for i in range(1, len(poly) - 1):
        tri = [poly[0], poly[i], poly[i + 1]]
        vals = [value(a, p) for p in tri]
        out += signed_area(tri) * (sum(v * v for v in vals) +
                sum(vals[j] * vals[k] for j, k in combinations(range(3), 2))) / 6
    return out


def convergence(n, L, eta, source=False):
    cells = exact_partition(n, L, eta, source)
    sup = F(0)
    grad_error = F(0)
    l2_error = F(0)
    grad_bound_squared = F(0)
    for poly, f, u in cells:
        diff = add(u, scale(-1, f))
        sup = max(sup, *(abs(value(diff, p)) for p in poly))
        grad_error += signed_area(poly) * (diff[0] ** 2 + diff[1] ** 2)
        l2_error += integral_square(poly, diff)
        grad_bound_squared = max(grad_bound_squared, u[0] ** 2 + u[1] ** 2)
        midpoint = tuple(sum(p[k] for p in poly) / len(poly) for k in range(2))
        for p in poly + [midpoint]:
            assert direct(n, L, F(0), p, source)[0] == value(f, p)
            assert direct(n, L, eta, p, source)[0] == value(u, p)
        assert direct(n, L, eta, midpoint, source)[1] == u[:2]
        assert direct(n, L, F(0), midpoint, source)[1] == f[:2]
    wl, wr, _ = weights(n, L, source)
    S = max(sum(wl), sum(wr))
    D = 2 * L * S
    assert sup <= 3 * D * eta
    assert grad_bound_squared <= 9 * (2 * S + 1) ** 2
    # Six graph bands and n wall strips: elementary area upper bound.
    exceptional_area_coefficient = 8 * n * L * L + 12 * D * (n + 1 + 2 * L)
    G = 3 * (2 * S + 1)
    area = 2 * L * (n + 1 + 2 * L)
    gradient_integral_coefficient = area * 36 * S * S + 4 * G * G * exceptional_area_coefficient
    assert grad_error <= gradient_integral_coefficient * eta
    assert l2_error <= area * sup ** 2
    return sup, grad_error, l2_error + grad_error, grad_bound_squared, len(cells)


def determinant(p, q, r):
    return (q[0] - p[0]) * (r[1] - p[1]) - (r[0] - p[0]) * (q[1] - p[1])


def incidence(n):
    b, delta, eps, _ = constants(n, 3)
    base = [(F(i), F(0) if i == 1 else b if i == 2 else eps * (n - i) ** 2)
            for i in range(1, n + 1)]
    points = [(x, y + j * delta) for x, y in base for j in range(1, 4)]
    ds = [abs(determinant(*triple)) for triple in combinations(base, 3)]
    assert min(ds) >= eps / (n + 1)
    assert all(d != 0 for d in ds)
    count = 0
    for triple in combinations(points, 3):
        if len({p[0] for p in triple}) != 1:
            assert determinant(*triple) != 0
            count += 1
    # Exact coverage of every possible nonvertical line through two designated points.
    lines = set()
    for p, q in combinations(points, 2):
        if p[0] != q[0]:
            slope = (q[1] - p[1]) / (q[0] - p[0])
            lines.add((slope, p[1] - slope * p[0]))
    assert max(sum(y == a * x + c for x, y in points) for a, c in lines) == 2
    # With l <= n lines, k designated walls, 3n <= 3k + 2(l-k) forces k=l=n.
    feasible = [(l, k) for l in range(n + 1) for k in range(l + 1)
                if 3 * n <= 3 * k + 2 * (l - k)]
    assert feasible == [(n, n)]
    wl, wr, _ = weights(n, 3)
    for i, k in base:
        hs = [max(value(h, (i, 0)), F(0)) for h in hidden(n, F(0))]
        side_weights = wl if i <= 2 else wr
        assert sum(a * h for a, h in zip(side_weights, hs)) == k
    print(f'incidence N={n}: {3*n} points; {count} nonvertical triples noncollinear; '
          f'minimum base determinant={min(ds)}; forced lines={feasible}')


def local_odd_part_grid(n):
    """Witness nonlinearity of the actual target odd part at all 3n bends."""
    b, delta, eps, lam = constants(n, 3)
    _, wr, _ = weights(n, 3)
    for i in range(1, n + 1):
        if i == 1:
            amin, aplus = F(0), b
        elif i == 2:
            amin, aplus = b, b + lam
        else:
            amin, aplus = -sum(wr[i-1:]), -sum(wr[i:])
        assert amin != aplus
        height = F(0) if i == 1 else b if i == 2 else eps * (n - i) ** 2
        t = delta / (100 * (1 + abs(amin) + abs(aplus)))
        ss = [min(amin, aplus) * t - t, min(amin, aplus) * t,
              (amin + aplus) * t / 2, max(amin, aplus) * t,
              max(amin, aplus) * t + t]
        for j in range(1, 4):
            P = (F(i), height + j * delta)
            odds = []
            for s in ss:
                plus = (P[0] + t, P[1] + s)
                minus = (P[0] - t, P[1] - s)
                odds.append(direct(n, 3, F(0), plus)[0] - direct(n, 3, F(0), minus)[0])
            slopes = {(odds[k+1] - odds[k]) / (ss[k+1] - ss[k]) for k in range(len(ss)-1)}
            assert len(slopes) > 1
    # A rational grid independently compares literal u_eta and the active-half expression.
    eta = F(1, 96)
    wl, wr, delta = weights(n, 3)
    for x, y in product([F(k, 4) for k in range(-12, 4 * (n + 4) + 1)],
                        [F(k, 4) for k in range(-12, 13)]):
        hs = [max(value(h, (x, y)), F(0)) for h in hidden(n, eta)]
        k = sum(a * h for a, h in zip(wl if x <= F(5, 2) else wr, hs))
        expected = sum(max(y - j * delta - k, F(0)) for j in range(1, 4))
        assert direct(n, 3, eta, (x, y))[0] == expected
    print(f'local odd-part and literal-network grid N={n}: OK')


def source_fidelity():
    # The source's exact four-wall family, not the asymmetric candidate family.
    # Reproduce K o S = h1 + lambda*h2 and Eq. (x2eta)'s active-half identity.
    n, L, eta = 4, 3, F(1, 96)
    lam = F(16)
    for x, y in product([F(k, 4) for k in range(-12, 33)], [F(k, 4) for k in range(-12, 13)]):
        hs = [max(value(h, (x, y)), F(0)) for h in hidden(n, eta)]
        mirrored = [max(value(h, (5 - x, y)), F(0)) for h in hidden(n, eta)]
        assert lam * mirrored[2] + mirrored[3] == hs[0] + lam * hs[1]
        if x >= F(5, 2):
            assert (hs[0] - hs[1] - 1) / eta == y
    b = [(F(1), F(0)), (F(2), F(1)), (F(3), F(1)), (F(4), F(0))]
    assert all(determinant(*t) for t in combinations(b, 3))
    for k in range(3):
        e = eta / (2 ** k)
        sup, grad2, w12, gb, cells = convergence(n, L, e, source=True)
        assert sup <= 3 * (lam + 2) * L * e
        print(f'source fidelity eta={e}: sup={sup}; gradient_L2_squared={grad2}; '
              f'W12_squared={w12}; cells={cells}')
    print('SOURCE_FIDELITY_OK: source equations, four-wall incidence, quantitative uniform bound')


def main():
    source_fidelity()
    for n in range(4, 13):
        incidence(n)
        local_odd_part_grid(n)
    previous = None
    for k in range(6):
        eta = F(1, 96 * 2 ** k)
        sup, grad2, w12, gb, cells = convergence(4, 3, eta)
        if previous is not None:
            assert sup < previous[0] and grad2 < previous[1] and w12 < previous[2]
        previous = (sup, grad2, w12)
        assert 41000 * eta < grad2 < 42000 * eta
        assert w12 < 43000 * eta
        assert gb < 2400
        print(f'candidate N=4 L=3 eta={eta}: sup={sup}; sup/eta={sup/eta}; '
              f'gradient_L2_squared={grad2}; gradient_L2_squared/eta={grad2/eta}; '
              f'W12_squared={w12}; max_gradient_squared={gb}; cells={cells}')
    print('N=4 six eta values: sup/eta=4617/32; 41000 < gradient_L2_squared/eta < 42000; W12_squared/eta < 43000; max_gradient_squared < 2400')
    print('Finite exact checks support the outlined general argument; they do not prove it.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
