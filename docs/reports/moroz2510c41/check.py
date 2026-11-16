#!/usr/bin/env python3
"""Exact arithmetic checks for Moroz's Perron continuity conjecture.

The checks use only integers and fractions.  For the constant-P specialization
P = 1, L is the positive Lüroth map on the branch containing x:
    L(x) = n(n+1)x - n,  1/(n+1) < x <= 1/n.
The formulas for the alternating coding are the affine maps
g_a(v) = 1/a - v/(a(a+1)).
"""

from fractions import Fraction


def L(x: Fraction) -> Fraction:
    """Positive Lüroth map, with its exact branch index."""
    for n in range(1, 10000):
        if Fraction(1, n + 1) < x <= Fraction(1, n):
            return n * (n + 1) * x - n
    raise AssertionError(f"no branch for {x}")


def g(a: int, v: Fraction) -> Fraction:
    return Fraction(1, a) - v / (a * (a + 1))


def r_family(a: int) -> Fraction:
    return Fraction(11 * a + 10, 11 * a * (a + 1))


def s_family(a: int) -> Fraction:
    return Fraction(11 * a + 1, 11 * a * (a + 1))


def main() -> None:
    # The source's explicit jump point: s_1 = 6/11, with F(s_1)=379/396.
    s1 = s_family(1)
    assert s1 == Fraction(6, 11)
    f_s1 = g(1, g(11, Fraction(2, 3)))
    assert f_s1 == Fraction(379, 396)

    # The stated right limit and its nonzero jump size.
    right_limit_s1 = Fraction(19, 20)
    assert right_limit_s1 - f_s1 == -Fraction(7, 990)
    assert right_limit_s1 == Fraction(10 * 1 + 9, 10 * 1 * (1 + 1))
    assert f_s1 - right_limit_s1 == Fraction(7, 990)

    # Verify both infinite families' exact formulas and the claimed nonzero gap.
    for a in range(1, 101):
        r = r_family(a)
        s = s_family(a)
        assert Fraction(0) < r <= Fraction(1)
        assert Fraction(0) < s <= Fraction(1)
        f_s = g(a, g(11, Fraction(2, 3)))
        right_s = Fraction(10 * a + 9, 10 * a * (a + 1))
        assert f_s == Fraction(198 * a + 181, 198 * a * (a + 1))
        assert right_s - f_s == -Fraction(7, 495 * a * (a + 1))
        # r_a has the positive Lüroth orbit stated in the settlement outline.
        orbit = [r]
        for _ in range(5):
            orbit.append(L(orbit[-1]))
        assert orbit[1:] == [
            Fraction(10, 11), Fraction(9, 11), Fraction(7, 11),
            Fraction(3, 11), Fraction(3, 11),
        ]

    print("s_1 =", s1)
    print("F(s_1) =", f_s1)
    print("right limit at s_1 =", right_limit_s1)
    print("verified r_a orbit and s_a jump formulas for a=1..100")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
