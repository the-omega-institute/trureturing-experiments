#!/usr/bin/env python3
"""Exact finite checks for Jenzer, arXiv:2610.21593v1, Remark 5.7 (TeX).

Sparse polynomials keep c_i and s_i independent. A monomial is encoded by
the exponents of s_i in base 3; the exponents of c_i are determined by the
number of occurrences of index i. No trigonometric evaluation or floats occur.
"""

from fractions import Fraction
from itertools import combinations, permutations, product


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def blade_mul(a, b, signature):
    """Geometric product in the orthogonal basis e_i^2 = signature[i]."""
    sign = 1
    for i, square in enumerate(signature):
        if (a >> i) & 1:
            if (b & ((1 << i) - 1)).bit_count() % 2:
                sign = -sign
            if (b >> i) & 1:
                sign *= square
    return a ^ b, sign


def tidy(d):
    return {key: value for key, value in d.items() if value}


def mv_mul(a, b, signature):
    out = {}
    for i, x in a.items():
        for j, y in b.items():
            blade, sign = blade_mul(i, j, signature)
            out[blade] = out.get(blade, 0) + x * y * sign
    return tidy(out)


def validate_family(family, signature):
    for i, b in enumerate(family):
        require(all(m.bit_count() == 2 for m in b), 'not a bivector')
        require(mv_mul(b, b, signature) == {0: -1}, 'wrong square')
        for a in family[:i]:
            ab = mv_mul(a, b, signature)
            ba = mv_mul(b, a, signature)
            require(ab == {m: -x for m, x in ba.items()}, 'not anticommuting')


def word_polynomial(family, word, signature):
    """Multiply (c_i + sign*s_i*B_i) factor by factor, in the given order."""
    polynomial = {(0, 0): 1}
    for i, factor_sign in word:
        out = {}
        for (blade, exponent), coefficient in polynomial.items():
            key = blade, exponent
            out[key] = out.get(key, 0) + coefficient
            for b, bc in family[i].items():
                new_blade, sign = blade_mul(blade, b, signature)
                key = new_blade, exponent + 3 ** i
                out[key] = out.get(key, 0) + coefficient * factor_sign * bc * sign
        polynomial = tidy(out)
    return polynomial


def scalar_coefficients(polynomial, blade=0):
    return {exponent: coefficient for (m, exponent), coefficient in polynomial.items() if m == blade}


def inversions(order, subset):
    positions = {j: i for i, j in enumerate(order)}
    return sum(positions[i] > positions[j] for i, j in combinations(subset, 2))


def expected_polynomial(k, order):
    """Independent subset enumeration of the source formula with c_i^2 weights."""
    out = {}
    for size in range(k + 1):
        for subset in combinations(range(k), size):
            exponent = sum(2 * 3 ** i for i in subset)
            out[exponent] = (-1) ** (size // 2 + inversions(order, subset))
    return out


def check_loop(family, signature, order):
    k = len(family)
    word = [(i, 1) for i in range(k)] + [(i, -1) for i in order]
    actual = scalar_coefficients(word_polynomial(family, word, signature))
    require(actual == expected_polynomial(k, order), ('scalar identity failed', k, order))


def fidelity():
    # Remark 4.4: B1=e1e2, B2=e1e3, B3=e2e3 in Cl_{0,3}.
    signature = (-1, -1, -1)
    triangle = [{3: 1}, {5: 1}, {6: 1}]
    validate_family(triangle, signature)
    p = mv_mul(mv_mul(triangle[0], triangle[1], signature), triangle[2], signature)
    require(p == {0: -1}, 'Remark 4.4 full product')
    chain = word_polynomial(triangle, [(i, 1) for i in range(3)], signature)
    require(scalar_coefficients(chain) == {0: 1, 13: -1}, 'Remark 4.4 chain polynomial')
    for order in permutations(range(3)):
        check_loop(triangle, signature, order)
    print('FIDELITY Remark 4.4: B1*B2*B3=-1; scalar(P3)=c1*c2*c3-s1*s2*s3')
    print('TRIANGLE all 6 return permutations: exact polynomial agreement; star condition fails')

    # An exact non-monomial bivector family with the same nonzero full trace.
    rotated = [{3: Fraction(3, 5), 5: Fraction(4, 5)},
               {3: Fraction(-4, 5), 5: Fraction(3, 5)}, {6: 1}]
    validate_family(rotated, signature)
    require(mv_mul(mv_mul(rotated[0], rotated[1], signature), rotated[2], signature)
            == {0: -1}, 'rotated full trace')
    for order in permutations(range(3)):
        check_loop(rotated, signature, order)
    print('RATIONAL ROTATION all 6 permutations: exact agreement for dependent non-monomial bivectors')


def universal_trace_checks():
    """All small permutations in the abstract Clifford presentation.

    This is stronger than using its canonical scalar projection alone: check
    the full-support coefficient separately for odd k, since a normalized
    tracial functional on a representation may retain this component.
    """
    total = 0
    for k in range(1, 7):
        signature = (-1,) * k
        formal_family = [{1 << i: 1} for i in range(k)]
        allowed = {0} | ({(1 << k) - 1} if k % 2 else set())
        # Every excluded monomial has an explicit conjugation sign witness.
        for d in range(1, 1 << k):
            if d in allowed:
                continue
            found = False
            for j in range(k):
                a, s1 = blade_mul(1 << j, d, signature)
                b, s2 = blade_mul(a, 1 << j, signature)
                require(b == d, 'conjugation support')
                if -s1 * s2 == -1:  # B_j^{-1}=-B_j
                    found = True
                    break
            require(found, ('missing trace-vanishing witness', k, d))
        count = 0
        for order in permutations(range(k)):
            word = [(i, 1) for i in range(k)] + [(i, -1) for i in order]
            actual = word_polynomial(formal_family, word, signature)
            require(scalar_coefficients(actual) == expected_polynomial(k, order),
                    ('diagonal coefficient', k, order))
            if k % 2:
                require(not scalar_coefficients(actual, (1 << k) - 1),
                        ('odd full-support cancellation', k, order))
            count += 1
        total += count
        print(f'UNIVERSAL k={k}: {count} permutations; independent c_i,s_i; all tracial coefficients agree')
    require(total == 873, 'permutation count')
    print(f'UNIVERSAL TOTAL {total} permutations for 1<=k<=6; odd exceptional coefficient exactly zero')


def geometric_small_cases():
    # Every signed orthogonal basis, every subset family of canonical bivectors,
    # and every order on that family, in all nondegenerate signatures for n<=4.
    family_count = permutation_count = 0
    for n in range(2, 5):
        blades = [(1 << i) | (1 << j) for i, j in combinations(range(n), 2)]
        for signature in product((-1, 1), repeat=n):
            units = [b for b in blades if blade_mul(b, b, signature) == (0, -1)]
            for k in range(1, len(units) + 1):
                for selection in combinations(units, k):
                    if not all(blade_mul(a, b, signature)[1] ==
                               -blade_mul(b, a, signature)[1]
                               for a, b in combinations(selection, 2)):
                        continue
                    family = [{b: 1} for b in selection]
                    validate_family(family, signature)
                    family_count += 1
                    for order in permutations(range(k)):
                        check_loop(family, signature, order)
                        permutation_count += 1
    print(f'GEOMETRIC n=2..4, all signatures and canonical bivector subset families: '
          f'{family_count} families, {permutation_count} orders, exact agreement')


if __name__ == '__main__':
    fidelity()
    universal_trace_checks()
    geometric_small_cases()
    print('LIMIT: finite exact checks support the preregistration; no unbounded theorem is certified here')
    print('ALL_CHECKS_OK')
