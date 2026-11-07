#!/usr/bin/env python3
"""Exact coefficient checks for arXiv:2608.30012v1, Section 3.

No transcendental numerical evaluations are used. Equality of restricted
linear-form multisets proves equality of their sinh products. At real test
points, nonzero (respectively positive) arguments certify nonzero sinh factors
using the mathematical facts sinh(t)=0 iff t=0 and sinh(t)>0 iff t>0.
These facts are not proved by this script; no Lean settlement is claimed.
"""

from collections import Counter
from fractions import Fraction
from itertools import combinations


# Table 1, equation coefficients (u,v,w), in printed row order.
VOGEL = (
    (1, 1, 0), (2, 1, 0), (1, 2, 0),
    (1, 0, 1), (0, 1, 1), (2, 0, 1), (0, 2, 1), (-2, -2, 1),
    (1, 0, 2), (0, 1, 2), (2, -1, 2), (-1, 2, 2),
)
A, B, C = 1000, 1000000, 1000000000


def subsets(items):
    for size in range(len(items) + 1):
        yield from combinations(items, size)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def restrict(form, normal):
    """Unique exact two-variable normal form on normal.dot(point)=0."""
    pivot = next(i for i in range(3) if normal[i])
    return tuple(Fraction(form[i]) - Fraction(form[pivot], normal[pivot]) * normal[i]
                 for i in range(3) if i != pivot)


def sinh_signature(forms, normal):
    """Account for sinh(-t)=-sinh(t), retaining multiplicities and zeros."""
    signature = Counter()
    sign = 1
    for f in forms:
        r = restrict(f, normal)
        nonzero = next((x for x in r if x), None)
        if nonzero is None:
            return 0, Counter()
        if nonzero < 0:
            sign = -sign
            r = tuple(-x for x in r)
        signature[r] += 1
    return sign, signature


def fidelity():
    # Equation (35), including the minus sign in -2 beta.
    red35 = ((0, -2, 0), (5, 2, 0), (4, 3, 0))
    green35 = ((1, 0, 0), (2, -1, 0), (6, 4, 0))
    for normal in VOGEL[:3]:
        assert sinh_signature(red35, normal) == sinh_signature(green35, normal)
        assert sinh_signature(green35, normal)[0] != 0
    print("FIDELITY equation (35): all three weight-zero line identities OK")

    # Equations (36)-(37): the source asserts equality of these multisets.
    red36 = ((1, 1, 0), (6, 5, 0), (5, 6, 0),
             (1, 0, 1), (0, 1, 1), (5, 5, 1))
    green37 = ((3, 2, 0), (2, 3, 0), (7, 7, 0),
               (-1, -1, 1), (4, 3, 1), (3, 4, 1))
    for normal in VOGEL[:3] + (VOGEL[7],):
        assert Counter(restrict(f, normal) for f in red36) == Counter(
            restrict(f, normal) for f in green37)
    print("FIDELITY equations (36)-(37): all four restricted multiset identities OK")


def construction(selected):
    normals = tuple(VOGEL[i] for i in selected) if selected else ((1, 0, 0),)
    forms = {}
    for mask in range(1 << len(normals)):
        chosen = [f for i, f in enumerate(normals) if mask & (1 << i)]
        a, b, c = (sum(f[j] for f in chosen) for j in range(3))
        forms[mask] = (A + a, B + b, c)
    red = [f for mask, f in forms.items() if mask.bit_count() % 2 == 0]
    green = [f for mask, f in forms.items() if mask.bit_count() % 2 == 1]
    return normals, forms, red, green


def check_converse():
    admissible = [s for s in subsets(tuple(range(12))) if sum(VOGEL[i][2] for i in s) < 3]
    assert len(admissible) == 160
    assert Counter(sum(VOGEL[i][2] for i in s) for s in admissible) == {0: 8, 1: 40, 2: 112}
    max_terms = 0
    min_denominator_argument = None
    line_identities = 0
    regular_restrictions = 0
    for selected in admissible:
        normals, forms, red, green = construction(selected)
        assert 1 <= len(normals) <= 5
        assert len(red) == len(green) == 1 << (len(normals) - 1)
        max_terms = max(max_terms, len(red))
        for a, b, c in forms.values():
            assert abs(a - A) <= 10 and abs(b - B) <= 10 and c in (0, 1, 2)

        # A direct parity-reversing bijection and, independently, normal forms.
        for j, index in enumerate(selected):
            normal = VOGEL[index]
            for mask, form in forms.items():
                other = mask ^ (1 << j)
                assert (mask.bit_count() - other.bit_count()) % 2 == 1
                assert restrict(form, normal) == restrict(forms[other], normal)
            assert Counter(restrict(f, normal) for f in red) == Counter(
                restrict(f, normal) for f in green)
            line_identities += 1

        # Stronger than the binding claim: both products have a nonzero real
        # value on every one of the twelve lines, including selected lines.
        for normal in VOGEL:
            u, v, w = normal
            point = (1, 0, -Fraction(u, w)) if w else (v, -u, 0)
            assert dot(normal, point) == 0
            for f in forms.values():
                assert any(restrict(f, normal))
                h = dot(f, point)
                assert h != 0
                if w:
                    assert h >= 986
                else:
                    assert h <= -997970
            regular_restrictions += 1

        witness = (B, -A, C)
        assert dot(forms[0], witness) == 0
        for f in green:
            h = dot(f, witness)
            assert h > 0
            min_denominator_argument = h if min_denominator_argument is None else min(
                min_denominator_argument, h)

    assert max_terms == 16
    assert min_denominator_argument == 998000
    print("CONVERSE admissible subsets: 160 (weights 0: 8; 1: 40; 2: 112)")
    print(f"CONVERSE selected-line parity/multiset identities: {line_identities}")
    print(f"CONVERSE nonzero numerator and denominator restrictions: {regular_restrictions}")
    print(f"CONVERSE maximum factors per side: {max_terms}")
    print(f"CONVERSE minimum denominator argument at nontriviality witness: {min_denominator_argument}")
    print("CONVERSE every numerator has a zero argument at the nontriviality witness")


def main():
    fidelity()
    check_converse()
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
