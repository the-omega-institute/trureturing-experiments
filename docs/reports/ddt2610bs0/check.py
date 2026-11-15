#!/usr/bin/env python3
"""Exact finite checks for the BS_0 nilring preregistration.

The carrier is the strictly upper triangular 3 by 3 matrices over F_2.  All
operations are integer arithmetic modulo 2; no floating point values occur.
"""

from itertools import product


MODULUS = 2
ZERO = ((0, 0, 0), (0, 0, 0), (0, 0, 0))
IDENTITY = ((1, 0, 0), (0, 1, 0), (0, 0, 1))


def add(a, b):
    return tuple(tuple((a[i][j] + b[i][j]) % MODULUS for j in range(3)) for i in range(3))


def mul(a, b):
    return tuple(
        tuple(sum(a[i][k] * b[k][j] for k in range(3)) % MODULUS for j in range(3))
        for i in range(3)
    )


def power(a, n):
    result = IDENTITY
    base = a
    while n:
        if n & 1:
            result = mul(result, base)
        base = mul(base, base)
        n >>= 1
    return result


def chain(xs):
    result = IDENTITY
    for x in xs:
        result = mul(result, x)
    return result


def strictly_upper_matrices():
    out = []
    for a01, a02, a12 in product(range(MODULUS), repeat=3):
        out.append(((0, a01, a02), (0, 0, a12), (0, 0, 0)))
    return out


R = strictly_upper_matrices()


def right_ann(a):
    return frozenset(i for i, x in enumerate(R) if mul(a, x) == ZERO)


def breadth_bound(max_family_size=4):
    principal = [right_ann(a) for a in R]
    for bound in range(1, len(R) + 1):
        for size in range(1, max_family_size + 1):
            for indices in product(range(len(R)), repeat=size):
                intersection = set(range(len(R)))
                for i in indices:
                    intersection &= principal[i]
                compressed = False
                for chosen_size in range(min(bound, size) + 1):
                    for chosen in product(indices, repeat=chosen_size):
                        candidate = set(range(len(R)))
                        for i in chosen:
                            candidate &= principal[i]
                        if candidate == intersection:
                            compressed = True
                            break
                    if compressed:
                        break
                if not compressed:
                    return None
        return bound


def second_order(a):
    return mul(a, a) == ZERO and all(mul(mul(a, x), a) == ZERO for x in R)


def main():
    assert len(R) == 8

    # Fidelity to Cherlin's Fact 3.14 for nilexponent n=3.
    assert all(power(x, 3) == ZERO for x in R)
    assert all(mul(mul(power(x, 2), y), power(x, 2)) == ZERO for x in R for y in R)

    # The source's Y_n family, with the formally adjoined identity in R^1.
    Y = {mul(b, power(x, 2)) for b in R + [IDENTITY] for x in R}
    assert all(second_order(y) for y in Y)
    for k in range(1, 8):
        if all(chain(xs) == ZERO for xs in product(sorted(Y, key=str), repeat=k)):
            y_product_length = k
            break
    else:
        raise AssertionError("no bounded product length for Y")

    # Directly check the expected settlement on this finite BS_0 example.
    breadth = breadth_bound()
    assert breadth is not None
    for nilpotence_class in range(1, 8):
        if all(chain(xs) == ZERO for xs in product(R, repeat=nilpotence_class)):
            break
    else:
        raise AssertionError("example was not nilpotent in the tested range")

    print("carrier_size=", len(R))
    print("verified_nilexponent=3")
    print("verified_second_order_product_length=", y_product_length)
    print("verified_right_BS0_breadth_bound=", breadth)
    print("verified_nilpotence_class=", nilpotence_class)
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
