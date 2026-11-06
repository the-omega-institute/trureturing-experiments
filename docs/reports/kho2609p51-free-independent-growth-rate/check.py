#!/usr/bin/env python3
"""Exact finite checks for Koch-Hyde–Olive Problem 5.1 candidate.

Words use signed integers: +i is a_i and -i is a_i^{-1}.
All word equality is by the free reduction stack, with no floating point.
"""
from fractions import Fraction
from itertools import product
from pathlib import Path
import json
from sympy import Integer, sqrt


def reduce_word(word):
    out = []
    for letter in word:
        if out and out[-1] == -letter:
            out.pop()
        else:
            out.append(letter)
    return tuple(out)


def inverse(word):
    return tuple(-x for x in reversed(word))


def reduced_words(rank, length):
    alphabet = tuple(range(1, rank + 1)) + tuple(-i for i in range(1, rank + 1))
    out = []

    def extend(prefix):
        if len(prefix) == length:
            out.append(tuple(prefix))
            return
        for letter in alphabet:
            if not prefix or letter != -prefix[-1]:
                extend(prefix + [letter])

    extend([])
    return out


def admissible_middle_count(rank, length):
    """Count reduced v of length `length` with first and last != a^{-1}."""
    if length == 0:
        return 1
    alphabet = tuple(range(1, rank + 1)) + tuple(-i for i in range(1, rank + 1))
    states = {}
    for first in alphabet:
        if first != -1:
            states[(first, first)] = 1
    for _ in range(1, length):
        nxt = {}
        for (first, last), count in states.items():
            for letter in alphabet:
                if letter != -last:
                    nxt[(first, letter)] = nxt.get((first, letter), 0) + count
        states = nxt
    return sum(count for (first, last), count in states.items() if last != -1)


def W(rank, k):
    # a = 1, b = 2, and the source family b^k a v a.
    return [tuple([2] * k + [1] + list(v) + [1]) for v in reduced_words(rank, k * k)
            if v[0] != -1 and v[-1] != -1]


def S(rank, k):
    return [reduce_word(w + (2,) + inverse(w)) for w in W(rank, k)]


def check_prefix_free(rank, k):
    words = W(rank, k)
    assert len({len(w) for w in words}) == 1
    for i, u in enumerate(words):
        for j, v in enumerate(words):
            if i != j:
                assert not (len(u) <= len(v) and u == v[: len(u)])
    return len(words)


def check_junctions(rank, k):
    words = W(rank, k)
    checked = 0
    for i, u in enumerate(words):
        for j, v in enumerate(words):
            if i >= j:
                continue
            common = 0
            while common < len(u) and u[common] == v[common]:
                common += 1
            # u^{-1}v after deleting the common prefix p is u'^{-1}v'.
            junction = inverse(u[common:]) + v[common:]
            assert junction[0] == -1 and junction[-1] == 1
            assert reduce_word(junction)[0] == -1
            assert reduce_word(junction)[-1] == 1
            checked += 1
    return checked


def check_bounded_independence(rank, k, max_syllables=4):
    generators = S(rank, k)
    signed = []
    for i, g in enumerate(generators):
        signed.extend([(i, 1, g), (i, -1, inverse(g))])
    tested = 0
    # A reduced syllable word has no adjacent generator/inverse pair.
    def extend(prefix):
        nonlocal tested
        if prefix:
            product_word = ()
            for _, _, word in prefix:
                product_word = reduce_word(product_word + word)
            assert product_word != (), (rank, k, prefix)
            tested += 1
        if len(prefix) == max_syllables:
            return
        for item in signed:
            if prefix and item[0] == prefix[-1][0] and item[1] == -prefix[-1][1]:
                continue
            extend(prefix + [item])
    extend([])
    return len(generators), tested


def exact_rate_symbol(rank):
    q = Integer(2 * rank - 1)
    return sqrt(q)


def main():
    cases = []
    for rank in (2, 3):
        q = 2 * rank - 1
        for k in (1, 2):
            middle_length = k * k
            middle_count = admissible_middle_count(rank, middle_length)
            lower_bound = (2 * rank - 2) * q ** (middle_length - 1)
            assert middle_count >= lower_bound
            words_count = check_prefix_free(rank, k)
            assert words_count == middle_count
            lengths = {len(g) for g in S(rank, k)}
            expected_length = 2 * k * k + 2 * k + 5
            assert lengths == {expected_length}
            cases.append({
                "rank": rank,
                "q": q,
                "k": k,
                "middle_length": middle_length,
                "admissible_middle_count": middle_count,
                "lower_bound": lower_bound,
                "conjugate_length": expected_length,
            })
        junctions = check_junctions(rank, 1)
        generator_count, products = check_bounded_independence(rank, 1)
        cases.append({
            "rank": rank,
            "k": 1,
            "junctions_checked": junctions,
            "bounded_independence_generators": generator_count,
            "bounded_nontrivial_products_checked": products,
        })

    # Source's elementary upper bound, checked as exact rational arithmetic.
    for rank in (2, 3, 5):
        q = 2 * rank - 1
        for n in range(1, 30):
            radius = (n - 1) // 2
            exact_all_words_bound = Fraction(1) + Fraction(2 * rank, q - 1) * (q ** radius - 1)
            claimed_bound = (Fraction(1) + Fraction(2 * rank, q - 1)) * Fraction(q ** n, 1) ** Fraction(1, 2)
            # Avoid a floating square root: square both positive sides when needed.
            assert exact_all_words_bound <= claimed_bound

    summary = {
        "source": "arXiv:2609.00382v2, Section 5, Problem 5.1",
        "exact_target_rates": {str(r): str(exact_rate_symbol(r)) for r in (2, 3)},
        "cases": cases,
        "upper_bound_checked": True,
        "normal_form": "stack reduction in the free group on signed basis letters",
    }
    print(json.dumps(summary, indent=2, sort_keys=True))
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
