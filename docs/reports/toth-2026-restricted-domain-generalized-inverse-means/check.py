#!/usr/bin/env python3
"""Exact arithmetic checks for Tóth, arXiv:2608.19997v1.

The checks use Fraction throughout.  They reproduce the empty-interior
obstruction (Example 2), verify an affine generator on a closed half-plane,
and exhibit the Jensen failure of a strictly increasing non-affine generator
on a closed disc.
"""

from fractions import Fraction as Q
from itertools import product


def dot(x, y):
    return sum(a * b for a, b in zip(x, y))


def add(x, y):
    return tuple(a + b for a, b in zip(x, y))


def scale(a, x):
    return tuple(a * z for z in x)


def midpoint(x, y):
    return scale(Q(1, 2), add(x, y))


def assert_true(condition, message):
    if not condition:
        raise AssertionError(message)


def check_example_2():
    # C=[0,1]×{0}, f(t,0)=(t,g(t)) with a deliberately non-affine g.
    g = lambda t: t**3 + t
    points = [(Q(i, 4), Q(0)) for i in range(5)]
    for x, y in product(points, repeat=2):
        fx = (x[0], g(x[0]))
        fy = (y[0], g(y[0]))
        if x != y:
            assert_true(dot(add(fx, scale(-1, fy)), add(x, scale(-1, y))) > 0,
                        "Example 2 strict increase")
        z = midpoint(x, y)
        fmean = scale(Q(1, 2), add(fx, fy))
        # The source's generalized inverse on conv(f(C)) is (s,0).
        F_fmean = (fmean[0], Q(0))
        assert_true(F_fmean == z, "Example 2 mean identity")
        assert_true((g(x[0]) != 0 or g(y[0]) != 0) or x == y,
                    "nonzero image witness for the ratio-free encoding")
    assert_true(g(Q(1, 2)) != (g(Q(0)) + g(Q(1))) / 2,
                "g is non-affine")


def check_affine_half_plane():
    # Closed half-plane C={ (x,y) : y>=0 }, A is positive definite as a
    # quadratic form but nonsymmetric, and b is nonzero.
    A = ((Q(1), Q(1)), (Q(0), Q(1)))
    b = (Q(2), Q(-1))

    def Af(v):
        return (A[0][0] * v[0] + A[0][1] * v[1],
                A[1][0] * v[0] + A[1][1] * v[1])

    def f(v):
        return add(Af(v), b)

    def Ainv(v):
        return (v[0] - v[1], v[1])

    def F(w):
        return Ainv(add(w, scale(-1, b)))

    # The quadratic form is x^2+xy+y^2=(x+y/2)^2+3y^2/4.
    for p, q in product(range(-3, 4), repeat=2):
        v = (Q(p), Q(q))
        if v != (Q(0), Q(0)):
            assert_true(dot(Af(v), v) > 0, "positive quadratic form")

    points = [(Q(i, 2), Q(j, 2)) for i in range(-2, 3) for j in range(5)]
    for xs in product(points, repeat=2):
        x, y = xs
        z = midpoint(x, y)
        lhs = F(midpoint(f(x), f(y)))
        assert_true(lhs == z, "affine half-plane mean identity")
        assert_true(f(z) == midpoint(f(x), f(y)), "affine Jensen identity")
    # Boundary propagation instance: x is interior, y is on the boundary.
    x, y, t = (Q(1), Q(2)), (Q(-1), Q(0)), Q(1, 3)
    z = add(scale(t, y), scale(1 - t, x))
    assert_true(f(z) == add(scale(t, f(y)), scale(1 - t, f(x))),
                "boundary propagation")


def check_non_affine_disc():
    # On the closed unit disc, f(u,v)=(u+u^3,v) is strictly increasing:
    # <f(x)-f(y),x-y> = (du)^2(1+u^2+uv+v^2)+(dv)^2 > 0.
    def f(p):
        u, v = p
        return (u + u**3, v)

    grid = [Q(i, 4) for i in range(-4, 5)]
    disc = [(u, v) for u, v in product(grid, repeat=2) if u * u + v * v <= 1]
    for x, y in product(disc, repeat=2):
        if x != y:
            du, dv = x[0] - y[0], x[1] - y[1]
            lhs = dot(add(f(x), scale(-1, f(y))), add(x, scale(-1, y)))
            factor = 1 + x[0] * x[0] + x[0] * y[0] + y[0] * y[0]
            assert_true(lhs == du * du * factor + dv * dv and lhs > 0,
                        "strict increase on disc grid")

    # The source conclusion would force Jensen affinity on interior points.
    x, y = (Q(0), Q(0)), (Q(1), Q(0))
    arithmetic_midpoint = midpoint(x, y)
    image_midpoint = midpoint(f(x), f(y))
    assert_true(f(arithmetic_midpoint) != image_midpoint,
                "non-affine disc violates the required Jensen identity")
    assert_true(f(arithmetic_midpoint) == (Q(5, 8), Q(0)), "exact value f(1/2,0)")
    assert_true(image_midpoint == (Q(1), Q(0)), "exact image midpoint")


def main():
    check_example_2()
    check_affine_half_plane()
    check_non_affine_disc()
    print("Example 2: exact strict increase and mean identity verified")
    print("Affine half-plane: exact positive-definite and boundary checks verified")
    print("Non-affine disc: exact strict increase and Jensen failure verified")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
