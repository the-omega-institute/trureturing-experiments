#!/usr/bin/env python3
"""Exact finite checks for Wołoszyn, arXiv:2605.14197v1, Question 5.5.

This program checks finite group witnesses, free-word identities, and the modal
cardinality pattern. It does not certify the infinitary HNN embedding theorem,
Cantor's theorem, existence of countable EC groups, or the universal refutation.
All mathematical computations use integers, tuples, sets and Boolean values.
"""
from itertools import permutations, product


def compose(p, q):
    return tuple(p[q[i]] for i in range(len(p)))


def cyclic(n):
    return tuple(tuple((i + a) % n for i in range(n)) for a in range(n))


def elementary_two(rank):
    return tuple(tuple(i ^ a for i in range(1 << rank)) for a in range(1 << rank))


def symmetric(n):
    return tuple(permutations(range(n)))


def closure(generators, identity):
    seen = {identity}
    frontier = [identity]
    while frontier:
        u = frontier.pop()
        for g in generators:
            v = compose(g, u)
            if v not in seen:
                seen.add(v)
                frontier.append(v)
    return seen


def group_axioms(elements):
    universe = set(elements)
    identity = tuple(range(len(elements[0])))
    assert identity in universe
    for a in elements:
        assert compose(identity, a) == a == compose(a, identity)
        assert any(compose(a, b) == identity == compose(b, a) for b in elements)
        for b in elements:
            assert compose(a, b) in universe
    # Associativity follows from actual function composition; check these tables too.
    assert all(compose(compose(a, b), c) == compose(a, compose(b, c))
               for a, b, c in product(elements, repeat=3))


def left_regular_double(elements):
    index = {g: i for i, g in enumerate(elements)}
    return {g: tuple(2 * index[compose(g, h)] + bit
                     for h in elements for bit in range(2)) for g in elements}


def membership_fidelity(elements):
    """Theorem 2.4: exact finite witnesses for every outside-subgroup element.

    Embed G in Sym(G x {0,1}) by doubled left translation. Swap the two
    copies precisely above A=<a,b>. This permutation centralizes A and fails
    to centralize every y outside A. This checks an actual overgroup, not
    merely the double centralizer inside the original group.
    """
    identity = tuple(range(len(elements[0])))
    regular = left_regular_double(elements)
    assert len(set(regular.values())) == len(elements)
    assert all(regular[compose(g, h)] == compose(regular[g], regular[h])
               for g, h in product(elements, repeat=2))
    witnesses = 0
    inside = 0
    for a, b in product(elements, repeat=2):
        subgroup = closure((a, b), identity)
        swap = tuple(2 * i + (1 - bit if h in subgroup else bit)
                     for i, h in enumerate(elements) for bit in range(2))
        assert set(swap) == set(range(2 * len(elements)))
        assert compose(swap, regular[a]) == compose(regular[a], swap)
        assert compose(swap, regular[b]) == compose(regular[b], swap)
        for y in elements:
            commutes = compose(swap, regular[y]) == compose(regular[y], swap)
            assert commutes == (y in subgroup)
            if y in subgroup:
                inside += 1
            else:
                witnesses += 1
    return inside, witnesses


def reduce_word(word):
    stack = []
    for letter in word:
        if stack and stack[-1] == -letter:
            stack.pop()
        else:
            stack.append(letter)
    return tuple(stack)


def inv(word):
    return tuple(-x for x in reversed(word))


def power(word, exponent):
    return reduce_word((word if exponent >= 0 else inv(word)) * abs(exponent))


def mul(*words):
    return reduce_word(x for word in words for x in word)


def free_basis_tests():
    # Project G*F(x,y) onto F(x,y). For i=0 these are x and y respectively.
    x, y = (1,), (2,)
    bases = [tuple(mul(power(y, -i), x, power(y, i)) for i in range(4)),
             tuple(mul(power(x, -i), y, power(x, i)) for i in range(4))]
    count = 0
    for basis in bases:
        frontier = [()]
        for length in range(1, 6):
            nxt = []
            for word in frontier:
                for k in range(1, 5):
                    for letter in (k, -k):
                        if word and word[-1] == -letter:
                            continue
                        new = word + (letter,)
                        image = mul(*(basis[abs(c)-1] if c > 0
                                      else inv(basis[abs(c)-1]) for c in new))
                        assert image, (length, new)
                        nxt.append(new)
                        count += 1
            frontier = nxt
    return count


def hnn_relation_tests():
    # t^-1 x t=y. The solved expression for g_i is a word in x,t.
    x, t = (1,), (3,)
    y = mul(inv(t), x, t)
    for i in range(1, 65):
        gi = mul(inv(x), power(y, i), t, power(x, -i), y,
                 power(x, i), inv(t), power(y, -i))
        assert all(abs(c) in (1, 3) for c in gi)
        lhs = mul(inv(t), power(y, -i), x, gi, power(y, i), t)
        rhs = mul(power(x, -i), y, power(x, i))
        assert lhs == rhs
    return 64


def modal_button_tests():
    # This is the quotient pattern for C=diamond Gen_2 after the mathematical
    # equivalence C <-> Countable. It is not a finite group counterexample.
    worlds = ('countable', 'uncountable')
    access = {'countable': worlds, 'uncountable': ('uncountable',)}
    C = {'countable': True, 'uncountable': False}
    U = {w: not C[w] for w in worlds}
    def diamond(p):
        return {w: any(p[v] for v in access[w]) for w in worlds}
    def box(p):
        return {w: all(p[v] for v in access[w]) for w in worlds}
    failure = {w: diamond(box(U))[w] and not U[w] for w in worlds}
    assert failure == {'countable': True, 'uncountable': False}
    assert box(U) == U
    # Cantor diagonal finite-prefix diagnostics, not a proof of uncountability.
    for n in range(1, 33):
        rows = [tuple(((i + 1) * (j + 3) + i // 2) % 2 for j in range(n))
                for i in range(n)]
        diagonal = tuple(1 - rows[i][i] for i in range(n))
        assert all(diagonal[i] != rows[i][i] for i in range(n))
    return failure


def main():
    groups = [('C1', cyclic(1)), ('C2', cyclic(2)), ('C3', cyclic(3)),
              ('C4', cyclic(4)), ('C5', cyclic(5)), ('C2^2', elementary_two(2)),
              ('C2^3', elementary_two(3)), ('S3', symmetric(3))]
    s3 = symmetric(3)
    assert any(compose(a, b) != compose(b, a) for a, b in product(s3, repeat=2))
    for name, group in groups:
        group_axioms(group)
        inside, outside = membership_fidelity(group)
        # Source §1: G embeds in G x S3 and the product is nonabelian.
        identity = tuple(range(len(group[0])))
        product_group = set(product(group, s3))
        assert len({(g, tuple(range(3))) for g in group}) == len(group)
        assert any((compose(identity, identity), compose(a, b)) !=
                   (compose(identity, identity), compose(b, a))
                   for a, b in product(s3, repeat=2))
        assert len(product_group) == 6 * len(group)
        print(f'FIDELITY {name}: |G|={len(group)}; membership_inside={inside}; '
              f'outside_extension_witnesses={outside}; nonabelian_product=OK')
    print('FREE_BASIS_PROJECTION reduced_words_checked=' + str(free_basis_tests()))
    print('HNN_TWO_GENERATOR solved_relations_checked=' + str(hnn_relation_tests()))
    print('MODAL_COUNTABILITY_PATTERN ' + str(modal_button_tests()))
    print('SCOPE: finite algebra and modal pattern checked; infinite bridges require proofs')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
