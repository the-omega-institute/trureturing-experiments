#!/usr/bin/env python3
"""Exact finite checks for Zenger's Section 4 embedding question.

Source: https://arxiv.org/abs/2605.00523v1, Sections 2 and 4.
Guard prior art: https://arxiv.org/abs/2603.05131v1, Definition 6.1.
These finite checks are evidence, not a proof of the unbounded assertion.
"""

from functools import lru_cache
from itertools import product
import json


BOT, P, Q = ("bot",), ("atom", 0), ("atom", 1)


def neg(p):
    return ("imp", p, BOT)


def em(p):
    return ("or", p, neg(p))


def subformulas(p):
    children = () if p[0] in ("bot", "atom") else (p[2],) if p[0] == "K" else p[1:]
    return {p}.union(*(subformulas(q) for q in children))


def size(p):
    if p[0] in ("bot", "atom"):
        return 1
    if p[0] == "K":
        return 1 + size(p[2])
    return 1 + sum(map(size, p[1:]))


def guard(p):
    parts = [em(q) for q in sorted(subformulas(p), key=repr)]
    result = parts[0]
    for q in parts[1:]:
        result = ("and", result, q)
    return ("C", result)


def translation(p):
    return ("imp", guard(p), p)


def relation(n, mask):
    return tuple(frozenset(v for v in range(n) if mask & (1 << (n * w + v))) for w in range(n))


def reflexive(r):
    return all(w in r[w] for w in range(len(r)))


def transitive(r):
    return all(r[v] <= r[w] for w in range(len(r)) for v in r[w])


def closure(r):
    result = [set(s) | {w} for w, s in enumerate(r)]
    for k in range(len(r)):
        for w in range(len(r)):
            if k in result[w]:
                result[w].update(result[k])
    return tuple(map(frozenset, result))


def evaluator(order, relations, valuation, domain):
    """Point semantics; restricting domain matters for both modal operators."""
    reach = closure(tuple(frozenset().union(*(r[w] for r in relations)) & domain
                          if w in domain else frozenset() for w in range(len(order))))

    @lru_cache(None)
    def truth(p, w):
        tag = p[0]
        if tag == "bot":
            return False
        if tag == "atom":
            return w in valuation[p[1]]
        if tag == "and":
            return truth(p[1], w) and truth(p[2], w)
        if tag == "or":
            return truth(p[1], w) or truth(p[2], w)
        if tag == "imp":
            return all(not truth(p[1], v) or truth(p[2], v) for v in order[w] & domain)
        if tag == "K":
            return all(truth(p[2], v) for v in relations[p[1]][w] & domain)
        if tag == "C":
            return all(truth(p[1], v) for v in reach[w] & domain)
        raise ValueError(p)

    return truth


def main():
    by_size = {1: {BOT, P, Q}}
    for k in (2, 3):
        f = set()
        for p in by_size[k - 1]:
            f.update({("K", 0, p), ("K", 1, p), ("C", p)})
        for i in range(1, k - 1):
            for p, q in product(by_size[i], by_size[k - 1 - i]):
                f.update((op, p, q) for op in ("and", "or", "imp"))
        by_size[k] = f
    formulas = set.union(*by_size.values())
    # Include nested implications, modal antecedents, and genuinely mixed paths.
    for p in (P, Q, ("K", 0, P), ("K", 1, Q), ("C", P)):
        formulas.update({em(p), ("C", em(p)), ("imp", neg(neg(p)), p),
                         ("imp", ("K", 0, p), ("C", p)),
                         ("C", ("imp", p, ("K", 1, p)))})
    formulas = sorted(formulas, key=repr)
    translated = {p: translation(p) for p in formulas}
    guards = {p: guard(p) for p in formulas}
    subs = {p: subformulas(p) for p in formulas}
    stats = {"formulas": len(formulas), "distinct_subformulas": len(set.union(*subs.values())),
             "models": {}, "root_checks": 0, "guarded_subformula_checks": 0,
             "source_monotonicity_checks": 0, "source_fixed_point_checks": 0,
             "classical_countermodel_checks": 0, "classical_equivalence_checks": 0}
    for p in formulas:
        assert size(translated[p]) <= 4 * (size(p) + 1) ** 2
    for n in (1, 2):
        worlds = frozenset(range(n))
        identity = tuple(frozenset({w}) for w in range(n))
        all_relations = [relation(n, m) for m in range(1 << (n * n))]
        orders = [o for o in all_relations if reflexive(o) and transitive(o)
                  and all(w == v or w not in o[v] for w in worlds for v in o[w])]
        counts = {"arbitrary": 0, "T": 0, "S4": 0}
        for order in orders:
            admissible = [r for r in all_relations
                          if all(r[v] <= r[w] for w in worlds for v in order[w])]
            upsets = [frozenset(w for w in worlds if m & (1 << w)) for m in range(1 << n)]
            upsets = [s for s in upsets if all(order[w] <= s for w in s)]
            for relations in product(admissible, repeat=2):
                is_t = all(map(reflexive, relations))
                is_s4 = is_t and all(map(transitive, relations))
                reach = closure(tuple(relations[0][w] | relations[1][w] for w in worlds))
                for valuation in product(upsets, repeat=2):
                    counts["arbitrary"] += 1
                    counts["T"] += int(is_t)
                    counts["S4"] += int(is_s4)
                    intuitionistic = evaluator(order, relations, valuation, worlds)
                    classical = evaluator(identity, relations, valuation, worlds)
                    restricted = {}
                    for y in worlds:
                        region = reach[y]
                        # Check the literal individual frame conditions after restriction.
                        if is_t:
                            assert all(w in r[w] & region for r in relations for w in region)
                        if is_s4:
                            assert all((r[v] & region) <= (r[w] & region)
                                       for r in relations for w in region for v in r[w] & region)
                        restricted[y] = evaluator(identity, relations, valuation, region)
                    for p in formulas:
                        for w in worlds:
                            # Source fidelity: Section 2 Monotonicity lemma.
                            for v in order[w]:
                                assert not intuitionistic(p, w) or intuitionistic(p, v)
                                stats["source_monotonicity_checks"] += 1
                            # Source fidelity: C p <-> p and everybody-knows C p.
                            cp = ("C", p)
                            assert intuitionistic(cp, w) == (intuitionistic(p, w) and all(
                                intuitionistic(cp, v) for r in relations for v in r[w]))
                            stats["source_fixed_point_checks"] += 1
                            assert classical(translated[p], w) == classical(p, w)
                            stats["classical_equivalence_checks"] += 1
                            if intuitionistic(guards[p], w):
                                for u in reach[w]:
                                    for q in subs[p]:
                                        assert intuitionistic(q, u) == restricted[w](q, u)
                                        stats["guarded_subformula_checks"] += 1
                            # All future witnesses of a failed translated implication.
                            assert intuitionistic(translated[p], w) == all(
                                not intuitionistic(guards[p], y) or restricted[y](p, y)
                                for y in order[w])
                            stats["root_checks"] += 1
                            if order == identity:
                                assert intuitionistic(guards[p], w)
                                assert intuitionistic(translated[p], w) == classical(p, w)
                                stats["classical_countermodel_checks"] += 1
        stats["models"][str(n)] = counts
    assert stats["models"]["2"] == {"arbitrary": 5554, "T": 328, "S4": 328}
    print(json.dumps(stats, indent=2, sort_keys=True))
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
