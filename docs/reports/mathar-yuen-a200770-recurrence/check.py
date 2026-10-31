#!/usr/bin/env python3
"""Exact checks for the corrected Mathar--Yuen conjecture in OEIS A200770.

Source: https://oeis.org/A200770/internal, revision 14.
The finite word checks are evidence for the proposed argument, not a Lean proof.
Run with Python 3 and SymPy 1.14.0: python3 check.py
"""
from collections import Counter
from fractions import Fraction
from itertools import product
from math import comb
import sympy as sp

SOURCE_TERMS = (
    3, 3, 2, 7, 9, 6, 22, 30, 20, 75, 105, 70, 266, 378, 252,
    966, 1386, 924, 3564, 5148, 3432, 13299, 19305, 12870, 50050,
    72930, 48620, 189618, 277134, 184756, 722228, 1058148, 705432,
    2762942, 4056234, 2704156, 10608612, 15600900, 10400600,
    40859500, 60174900,
)
SOURCE_N4 = {
    (0, 0, 1, 2), (0, 2, 1, 2), (0, 0, 2, 1), (0, 1, 2, 2),
    (0, 1, 2, 1), (0, 1, 1, 2), (0, 2, 1, 1),
}


def require(condition, context):
    if not condition:
        raise AssertionError(context)


def source_averages(word):
    """Literal nondecreasing nonempty prefix averages, using rationals."""
    prefix = 0
    previous = None
    for length, letter in enumerate(word, start=1):
        prefix += letter
        average = Fraction(prefix, length)
        if previous is not None and average < previous:
            return False
        previous = average
    return True


def content(word):
    return tuple(word.count(letter) for letter in range(3))


def balanced(counts):
    return all(abs(r - s) <= 1 for r in counts for s in counts)


def classification(word):
    r, s, t = content(word)
    middle_twos = min(r, t)
    middle = word[r:r + s + middle_twos]
    return (
        word[:r] == (0,) * r
        and middle.count(1) == s
        and middle.count(2) == middle_twos
        and word[r + s + middle_twos:] == (2,) * (t - middle_twos)
    )


def fixed_content(r, s, t):
    return comb(s + min(r, t), s)


def residue_formula(n):
    require(n >= 1, ('offset', n))
    q, r = divmod(n, 3)
    c, b = comb(2 * q, q), comb(2 * q + 1, q)
    return (c, 2 * c + b, 3 * b)[r]


def recurrence_terms(n, a):
    return (
        (n + 2) * (2770 * n - 4893) * a(n),
        3 * (1090 * n**2 - 4919 * n + 1654) * a(n - 1),
        4686 * a(n - 2),
        2 * (-5540 * n**2 + 7016 * n + 17127) * a(n - 3),
        -6 * (1090 * n - 2739) * (2 * n - 7) * a(n - 4),
    )


def check_words():
    brute = {}
    tested = 0
    triples = 0
    for n in range(12):
        actual_counts = Counter()
        accepted = set()
        for word in product(range(3), repeat=n):
            actual = source_averages(word)
            require(actual == classification(word), ('classification', word))
            if actual:
                actual_counts[content(word)] += 1
                if balanced(content(word)):
                    accepted.add(word)
            tested += 1
        for r in range(n + 1):
            for s in range(n - r + 1):
                t = n - r - s
                require(actual_counts[r, s, t] == fixed_content(r, s, t),
                        ('fixed-content', r, s, t, actual_counts[r, s, t]))
                triples += 1
        if n:
            brute[n] = len(accepted)
            require(brute[n] == SOURCE_TERMS[n - 1], ('source data', n))
            require(brute[n] == residue_formula(n), ('residue formula', n))
        if n == 4:
            require(accepted == SOURCE_N4, ('source n=4 examples', accepted))
    for n in range(5, 12):
        require(sum(recurrence_terms(n, brute.__getitem__)) == 0,
                ('brute recurrence', n))
    print('EXHAUSTIVE_WORDS n=0..11:', tested)
    print('FIXED_CONTENT_TRIPLES total_length=0..11:', triples)
    print('SOURCE_N4: all seven arrays reproduced exactly')
    print('BRUTE_SOURCE_TERMS n=1..11:', ','.join(map(str, brute.values())))
    print('BRUTE_RECURRENCE n=5..11: all residuals zero')


def check_source_and_large_indices():
    for n, expected in enumerate(SOURCE_TERMS, 1):
        require(residue_formula(n) == expected, ('source formula', n))
    for n in range(1, 501):
        fixed_sum = sum(fixed_content(r, s, n - r - s)
                        for r in range(n // 3, n // 3 + 2)
                        for s in range(n // 3, n // 3 + 2)
                        if n - r - s >= 0 and balanced((r, s, n - r - s)))
        require(fixed_sum == residue_formula(n), ('balanced sum', n))
    for n in range(5, len(SOURCE_TERMS) + 1):
        require(sum(recurrence_terms(n, lambda k: SOURCE_TERMS[k - 1])) == 0,
                ('source recurrence', n))
    for n in range(5, 501):
        require(sum(recurrence_terms(n, residue_formula)) == 0,
                ('large recurrence', n))
    print('SOURCE_DATA n=1..41: all 41 printed terms matched')
    print('BALANCED_CONTENT_SUM n=1..500: agrees with three residue formulas')
    print('SOURCE_RECURRENCE n=5..41: all residuals zero')
    print('FORMULA_RECURRENCE n=5..500: all residuals zero')


def check_symbolic():
    q = sp.Symbol('q', integer=True, positive=True)
    central_ratio = {
        0: sp.Integer(1),
        -1: q / (2 * (2 * q - 1)),
        -2: q * (q - 1) / (4 * (2 * q - 1) * (2 * q - 3)),
    }

    def ratio(residue, shift):
        delta, r = divmod(residue - shift, 3)
        m = q + delta
        factors = (sp.Integer(1), (4 * m + 3) / (m + 1),
                   3 * (2 * m + 1) / (m + 1))
        return sp.cancel(factors[r] * central_ratio[delta])

    for residue, qmin in [(0, 2), (1, 2), (2, 1)]:
        n = 3 * q + residue
        coefficients = (
            (n + 2) * (2770 * n - 4893),
            3 * (1090 * n**2 - 4919 * n + 1654),
            sp.Integer(4686),
            2 * (-5540 * n**2 + 7016 * n + 17127),
            -6 * (1090 * n - 2739) * (2 * n - 7),
        )
        denominator = (4 * (2 * q - 1) if residue == 0
                       else 2 * (q + 1) * (2 * q - 1))
        polynomials = []
        print('SYMBOLIC_RESIDUE r=' + str(residue) + ' q>=' + str(qmin)
              + ' denominator=' + str(denominator))
        for shift, coefficient in enumerate(coefficients):
            term = sp.cancel(denominator * coefficient * ratio(residue, shift))
            numerator, divisor = sp.fraction(term)
            require(divisor == 1, ('polynomial denominator', residue, shift))
            poly = sp.Poly(sp.expand(numerator), q, domain=sp.ZZ)
            polynomials.append(poly)
            print('  P' + str(shift) + '=' + str(poly.as_expr()))
        total = sum(polynomials, sp.Poly(0, q, domain=sp.ZZ))
        require(total.is_zero, ('symbolic residue', residue, total))
        degree = max(poly.degree() for poly in polynomials)
        sums = [sum(poly.nth(k) for poly in polynomials)
                for k in range(degree, -1, -1)]
        print('  coefficient_sums_descending=' + str(sums))
        print('  residual=' + str(total.as_expr()))
    print('SYMBOLIC_CERTIFICATES: 15 integer polynomials, 3 zero identities')


def main():
    print('SOURCE: OEIS A200770 revision 14, corrected coefficient 2770')
    print('ARITHMETIC: Python integers, Fraction, SymPy exact ZZ polynomials')
    check_words()
    check_source_and_large_indices()
    check_symbolic()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
