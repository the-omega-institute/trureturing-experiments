#!/usr/bin/env python3
"""Exact finite checks for the Liu--Long--Shi--Su--Xiao Smith-set question.

All probabilities are Fraction values; no floating point is used.
"""
from fractions import Fraction
from itertools import permutations, product
from math import factorial


def luce_probability(perm, weights):
    remaining = list(perm)
    out = Fraction(1, 1)
    for item in perm:
        denom = sum(weights[j] for j in remaining)
        out *= weights[item] / denom
        remaining.remove(item)
    return out


def ranking_table(n, light=Fraction(1, 6)):
    weights = [light] + [Fraction(1, 1)] * (n - 1)
    rows = []
    for perm in permutations(range(n)):
        rows.append((perm, luce_probability(perm, weights)))
    assert sum(p for _, p in rows) == 1
    assert all(p > 0 for _, p in rows)
    return rows


def beats_mask(perm, n):
    pos = [0] * n
    for k, v in enumerate(perm):
        pos[v] = k
    mask = [0] * n
    for i in range(n):
        for j in range(n):
            if i != j and pos[i] < pos[j]:
                mask[i] |= 1 << j
    return mask


def smith_size(masks, n):
    # The Smith set is the smallest nonempty set that beats every vertex outside.
    for size in range(1, n + 1):
        for subset in range(1, 1 << n):
            if subset.bit_count() != size:
                continue
            outside = ((1 << n) - 1) ^ subset
            if all((masks[v] & outside) == outside for v in range(n) if subset & (1 << v)):
                return size
    raise AssertionError("nonempty dominating set must exist")


def exact_model_check(n):
    rows = ranking_table(n)
    perms = [p for p, _ in rows]
    probs = [p for _, p in rows]
    masks = [beats_mask(p, n) for p in perms]

    # Fidelity: under the Luce formula the light alternative beats a core item with 1/7.
    p_light_beats_one = sum(
        prob for perm, prob in zip(perms, probs) if perm.index(0) < perm.index(1)
    )
    assert p_light_beats_one == Fraction(1, 7)

    # Expected refutation witness: the light alternative is a Condorcet loser.
    p_loser = Fraction(0, 1)
    p_size_n_minus_1 = Fraction(0, 1)
    p_full_failure = Fraction(0, 1)
    for choices in product(range(len(perms)), repeat=3):
        weight = probs[choices[0]] * probs[choices[1]] * probs[choices[2]]
        all_masks = [masks[i] for i in choices]
        majority = [0] * n
        for i in range(n):
            for j in range(n):
                if i != j and sum(1 for mask in all_masks if mask[i] & (1 << j)) >= 2:
                    majority[i] |= 1 << j
        light_loses = all(not (majority[0] & (1 << j)) for j in range(1, n))
        if light_loses:
            p_loser += weight
        size = smith_size(majority, n)
        if size == n - 1:
            p_size_n_minus_1 += weight
        if size < n:
            p_full_failure += weight

    # The proposed proof's lower bound is (1/2)(6n)^(-1/4). Compare fourth powers exactly.
    lower_bound_fourth = Fraction(1, 96 * n)
    assert p_loser ** 4 >= lower_bound_fourth
    assert p_full_failure >= p_loser
    assert p_full_failure > 0
    assert p_size_n_minus_1 > 0
    return {
        "n": n,
        "ranking_count": factorial(n),
        "p_light_beats_core": p_light_beats_one,
        "p_condorcet_loser_0": p_loser,
        "p_smith_n_minus_1": p_size_n_minus_1,
        "p_full_failure": p_full_failure,
        "lower_bound_fourth": lower_bound_fourth,
    }


def main():
    results = [exact_model_check(n) for n in (4,)]
    for row in results:
        print("n={n} rankings={ranking_count}".format(**row))
        print("P(light beats a specified core item) = {p_light_beats_core}".format(**row))
        print("P(light is Condorcet loser) = {p_condorcet_loser_0}".format(**row))
        print("P(|S|=n-1) = {p_smith_n_minus_1}".format(**row))
        print("P(|S|<n) = {p_full_failure}".format(**row))
        print("exact lower-bound check: P(loser)^4 >= {lower_bound_fourth}".format(**row))
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
