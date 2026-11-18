#!/usr/bin/env python3
"""Exact fidelity and counterexample checks for arXiv:2603.23614v2."""
from collections import Counter
from fractions import Fraction
from itertools import permutations, product
from math import prod


def compose(left, right):
    return tuple(left[i] for i in right)


def inverse(perm):
    result = [0] * len(perm)
    for i, image in enumerate(perm):
        result[image] = i
    return tuple(result)


def even(perm):
    inversions = sum(a > b for i, a in enumerate(perm) for b in perm[i + 1:])
    return inversions % 2 == 0


def commutator(left, right):
    return compose(compose(left, right), inverse(compose(right, left)))


def commutator_counts(group):
    return Counter(commutator(x, y) for x in group for y in group)


def check_a5():
    identity = tuple(range(5))
    group = [g for g in permutations(identity) if even(g)]
    assert len(group) == 60
    sylow_five = set()
    for g in group:
        powers = [identity]
        for _ in range(4):
            powers.append(compose(powers[-1], g))
        if len(set(powers)) == 5 and compose(powers[-1], g) == identity:
            sylow_five.add(frozenset(powers))
    assert len(sylow_five) == 6

    # Remove the constant subrepresentation of the action on Sylow 5-subgroups.
    character = {}
    for g in group:
        g_inv = inverse(g)
        character[g] = sum(
            frozenset(compose(compose(g, h), g_inv) for h in subgroup) == subgroup
            for subgroup in sylow_five
        ) - 1
    assert character[identity] == 5
    assert sum(value**2 for value in character.values()) == len(group)
    assert Counter(map(abs, character.values())) == {5: 1, 1: 35, 0: 24}

    counts = commutator_counts(group)
    probabilities = Counter()
    for h, count in counts.items():
        probabilities[abs(character[h])] += Fraction(count, len(group)**2)
    assert probabilities == {5: Fraction(5, 60), 1: Fraction(29, 60), 0: Fraction(26, 60)}
    print("A5 Example 5.1: level-set sizes (1, 35, 24), probabilities (5, 29, 26)/60")

    for g in group:
        values = {h: abs(character[compose(inverse(h), g)]) for h in counts}
        average = Fraction(sum(counts[h] * value for h, value in values.items()), len(group)**2)
        gap = average / 5 - Fraction(counts[g], len(group)**2)
        tail = Fraction(sum(counts[h] * value for h, value in values.items() if h != g),
                        5 * len(group)**2)
        assert gap == tail
        assert sum(counts[h] for h in counts if values[h] == 5) == counts[g]
    print("A5 Lemmas 5.3 and 5.5(i) at t=1: checked for all 60 targets")


def check_s4_powers():
    identity = tuple(range(4))
    group = list(permutations(identity))
    character = {g: sum(i == g[i] for i in identity) - 1 for g in group}
    counts = commutator_counts(group)
    commutators = set(counts)
    assert len(group) == 24
    assert counts.total() == len(group)**2
    assert commutators == {g for g in group if even(g)}
    assert len(commutators) == 12
    assert Counter(character.values()) == {3: 1, 1: 6, -1: 9, 0: 8}
    assert sum(value**2 for value in character.values()) == len(group)
    assert [g for g in group if all(compose(g, h) == compose(h, g) for h in group)] == [identity]
    transposition = (1, 0, 2, 3)
    double_transposition = (1, 0, 3, 2)
    assert commutator((0, 2, 1, 3), transposition) == (1, 2, 0, 3)
    assert commutator((1, 2, 3, 0), (2, 1, 0, 3)) == (2, 3, 0, 1)

    # Summing over independent source pairs factors the unnormalized expectation.
    moments = {}
    for g in group:
        moments[g] = sum(count * abs(character[compose(inverse(h), g)])
                         for h, count in counts.items())
        direct = sum(abs(character[compose(inverse(commutator(x, y)), g)])
                     for x in group for y in group)
        assert direct == moments[g]

    for k in range(1, 4):
        degree = 3**k
        denominator = degree * len(group)**(2 * k)
        identity_power = (identity,) * k
        minimum_gap = None
        moduli = set()
        square_sum = 0
        supported_targets = 0
        for g in product(group, repeat=k):
            modulus = abs(prod(character[gi] for gi in g))
            moduli.add(modulus)
            square_sum += modulus**2
            assert (modulus == degree) == (g == identity_power)

            multiplicity = prod(counts[gi] for gi in g)
            supported_targets += multiplicity > 0
            gap = Fraction(prod(moments[gi] for gi in g) - degree * multiplicity, denominator)

            t = tuple(double_transposition if even(gi) else transposition for gi in g)
            h = tuple(compose(gi, inverse(ti)) for gi, ti in zip(g, t))
            assert h != g and all(hi in commutators for hi in h)
            assert tuple(compose(inverse(hi), gi) for hi, gi in zip(h, g)) == t
            assert abs(prod(character[ti] for ti in t)) == 1
            assert gap >= Fraction(1, denominator)
            minimum_gap = gap if minimum_gap is None else min(minimum_gap, gap)

        assert denominator == 1728**k
        assert square_sum == len(group)**k
        assert sorted(moduli)[-2] == 3**(k - 1)
        assert supported_targets == 12**k < 24**k
        print(f"S4^{k}: degree={degree}, gamma={3**(k - 1)}, "
              f"commutators={supported_targets}/{24**k}, min gap={minimum_gap}")


if __name__ == "__main__":
    check_a5()
    check_s4_powers()
    print("ALL_CHECKS_OK")
