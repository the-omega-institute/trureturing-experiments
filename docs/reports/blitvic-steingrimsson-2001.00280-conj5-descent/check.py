#!/usr/bin/env python3
"""Exact checks for the unresolved conventions of arXiv:2001.00280v4, Conjecture 5.

The internal-comparison variant is tested separately from Definition 5. No
settlement of the published conjecture is asserted. All evidence is integral.
"""
from collections import Counter
from functools import lru_cache
from itertools import combinations, permutations, product
from math import comb, factorial


def arrangements(n, k):
    for perm in permutations(range(1, n + 1)):
        fixed = [i for i, value in enumerate(perm, 1) if i == value]
        for colors in product(range(1, k + 1), repeat=len(fixed)):
            yield perm, dict(zip(fixed, colors))


def forms(n, k, derangement=False):
    alphabet = tuple(range(-k if derangement else 1 - k, 0))
    for negative_count in range(n + 1):
        m = n - negative_count
        for positive in permutations(range(1, m + 1)):
            if derangement and any(i == x for i, x in enumerate(positive, 1)):
                continue
            for positions in combinations(range(n), negative_count):
                for negative in product(alphabet, repeat=negative_count):
                    assigned = dict(zip(positions, negative))
                    pos = iter(positive)
                    yield tuple(assigned[i] if i in assigned else next(pos) for i in range(n))


def word_form(perm, colors, k, derangement=False):
    raw = [-colors[i] if i in colors and (derangement or colors[i] != k) else x
           for i, x in enumerate(perm, 1)]
    return positive_reduce(raw)


def positive_reduce(word):
    ranks = {x: i for i, x in enumerate(sorted(x for x in word if x > 0), 1)}
    return tuple(ranks[x] if x > 0 else x for x in word)


def descents(word):
    return sum(x > y for x, y in zip(word, word[1:]))


def colored_descents(perm, colors, k, convention, terminal=False):
    def color(i):
        if i not in colors:
            return 0
        if convention == 'literal':
            return colors[i]
        if convention == 'zero_based':
            return colors[i] - 1
        assert convention == 'cyclic' and k > 0
        return colors[i] % k
    pairs = tuple((color(i), x) for i, x in enumerate(perm, 1))
    if terminal:
        pairs += ((0, len(perm) + 1),)
    return descents(pairs)


def internal_bijection(perm, colors, k, convention):
    """Reverse-complement the permutation, mirror colors, then positive-reduce.

    Positive colors all become negative letters. With literal labels the image
    is a derangement form; with either zero-based labeling it is a permutation
    form. A zero-colored fixed point remains a positive letter.
    """
    n = len(perm)
    pairs = []
    for i in range(n, 0, -1):
        c = colors.get(i, 0)
        if i in colors and convention == 'zero_based':
            c -= 1
        if i in colors and convention == 'cyclic':
            c %= k
        pairs.append(-c if c > 0 else n + 1 - perm[i - 1])
    return positive_reduce(pairs)


def coefficients(hist, n):
    return tuple(hist[i] for i in range(n + 1))


def signatures(n):
    """Enumerate S_n exactly, retaining its fixed-point mask and descent mask.

    These masks determine every colored adjacent comparison, because equal
    colors use the original value comparison and unequal colors use colors.
    """
    out = Counter()
    for p in permutations(range(1, n + 1)):
        fixed = sum(1 << i for i, x in enumerate(p) if x == i + 1)
        desc = sum(1 << i for i in range(n - 1) if p[i] > p[i + 1])
        out[fixed, desc] += 1
    assert sum(out.values()) == factorial(n)
    return out


def color_transfer(n, fixed, desc, labels):
    """Sum all fixed-point colorings of one signature by an integer transfer."""
    if n == 0:
        return (1,), (1,)
    choices = labels if fixed & 1 else (0,)
    state = {c: [1] + [0] * n for c in choices}
    for i in range(1, n):
        choices = labels if fixed & (1 << i) else (0,)
        new = {c: [0] * (n + 1) for c in choices}
        base_descent = bool(desc & (1 << (i - 1)))
        for left, poly in state.items():
            for right in choices:
                shift = int(left > right or (left == right and base_descent))
                target = new[right]
                for d, count in enumerate(poly[:n + 1 - shift]):
                    target[d + shift] += count
        state = new
    internal, terminal = [0] * (n + 1), [0] * (n + 1)
    for last_color, poly in state.items():
        shift = int(last_color > 0)
        for d, count in enumerate(poly):
            internal[d] += count
            if count:
                terminal[d + shift] += count
    return tuple(internal), tuple(terminal)


def colored_polynomials(n, k, sigs, zero_based):
    labels = tuple(range(k)) if zero_based else tuple(range(1, k + 1))
    internal, terminal = [0] * (n + 1), [0] * (n + 1)
    for (fixed, desc), multiplicity in sigs.items():
        a, b = color_transfer(n, fixed, desc, labels)
        for d in range(n + 1):
            internal[d] += multiplicity * a[d]
            terminal[d] += multiplicity * b[d]
    return tuple(internal), tuple(terminal)


def word_polynomial(n, k, derangement):
    """Independent word enumeration, memoized by unused positive letters.

    Negative letters may repeat; each positive letter is used once. In a
    derangement form its j-th positive letter is forbidden to equal j.
    """
    r = k if derangement else k - 1
    assert r >= 0
    total = [0] * (n + 1)
    for m in range(n + 1):
        negative_total = n - m

        @lru_cache(None)
        def visit(unused, neg_left, last):
            if unused == 0 and neg_left == 0:
                return (1,) + (0,) * n
            out = [0] * (n + 1)
            choices = []
            if neg_left:
                choices.extend((x, unused, neg_left - 1) for x in range(-r, 0))
            pos_index = m - unused.bit_count() + 1
            for j in range(m):
                if unused & (1 << j) and (not derangement or j + 1 != pos_index):
                    choices.append((j + 1, unused ^ (1 << j), neg_left))
            for x, remaining, negatives in choices:
                shift = int(last != 0 and last > x)
                tail = visit(remaining, negatives, x)
                for d, count in enumerate(tail[:n + 1 - shift]):
                    out[d + shift] += count
            return tuple(out)

        result = visit((1 << m) - 1, negative_total, 0)
        for d, count in enumerate(result):
            total[d] += count
    return tuple(total)


def main():
    # Fidelity: source Proposition 2(1) and its permanent interpretation (3).
    recurrence = {}
    for k in range(5):
        for n in range(10):
            recurrence[n, k] = 1 if n == 0 else n * recurrence[n - 1, k] + (k - 1) ** n
            permanent = sum(k ** sum(i == x for i, x in enumerate(p, 1))
                            for p in permutations(range(1, n + 1)))
            assert permanent == recurrence[n, k]
        print('FIDELITY Proposition 2 k=', k, 'A_k(0..9)=',
              tuple(recurrence[n, k] for n in range(10)), flush=True)

    example = (6, 2, 1, 4, 5, 7, 3)
    coloring = {2: 4, 4: 1, 5: 4}
    assert word_form(example, coloring, 4, True) == (3, -4, 1, -1, -4, 4, 2)
    assert word_form(example, coloring, 4) == (5, 2, 1, -1, 4, 6, 3)
    print('FIDELITY source 4-arrangement example OK', flush=True)

    # The type-B row in the source table has alpha_0 = u = 1+t.
    signed_one = [((1,), {}), ((1,), {1: 2})]
    b1_terminal = Counter(colored_descents(p, c, 2, 'zero_based', True)
                          for p, c in signed_one)
    b1_internal = Counter(colored_descents(p, c, 2, 'zero_based')
                          for p, c in signed_one)
    assert coefficients(b1_terminal, 1) == (1, 1)
    assert coefficients(b1_internal, 1) == (2, 0)
    print('FIDELITY type-B size one: terminal=(1,1), internal=(2,0)', flush=True)

    # Explicit sets independently cross-check both DPs and the candidate maps.
    for n in range(6):
        sigs = signatures(n)
        for k in range(1, 5):
            objects = tuple(arrangements(n, k))
            pwords, dwords = set(forms(n, k)), set(forms(n, k, True))
            assert len(pwords) == len(dwords) == len(objects) == recurrence[n, k]
            assert {word_form(p, c, k) for p, c in objects} == pwords
            assert {word_form(p, c, k, True) for p, c in objects} == dwords
            assert coefficients(Counter(map(descents, pwords)), n) == word_polynomial(n, k, False)
            assert coefficients(Counter(map(descents, dwords)), n) == word_polynomial(n, k, True)
            for convention in ('literal', 'zero_based', 'cyclic'):
                image = set()
                hist, terminal = Counter(), Counter()
                for p, c in objects:
                    w = internal_bijection(p, c, k, convention)
                    d = colored_descents(p, c, k, convention)
                    assert descents(w) == d
                    image.add(w)
                    hist[d] += 1
                    terminal[colored_descents(p, c, k, convention, True)] += 1
                assert image == (dwords if convention == 'literal' else pwords)
                dp, tp = colored_polynomials(n, k, sigs, convention != 'literal')
                assert coefficients(hist, n) == dp
                assert coefficients(terminal, n) == tp
    print('EXPLICIT_CROSSCHECK 0<=n<=5, 1<=k<=4: sets, transfer DPs, bijections OK', flush=True)

    # Full finite distributions, without sampling or floating-point arithmetic.
    for n in range(10):
        sigs = signatures(n)
        print('SIGNATURES n=', n, 'permutations=', factorial(n), 'states=', len(sigs), flush=True)
        for k in range(1, 5):
            positive, literal_terminal = colored_polynomials(n, k, sigs, False)
            shifted, shifted_terminal = colored_polynomials(n, k, sigs, True)
            # phi-1 and phi mod k each biject [k] onto the same color alphabet.
            pword = word_polynomial(n, k, False)
            dword = word_polynomial(n, k, True)
            assert positive == shifted == pword == dword
            assert sum(pword) == recurrence[n, k]
            assert pword[0] == comb(n + k - 1, k - 1)
            if n:
                assert literal_terminal[0] == 0 and shifted_terminal[0] == 1
            if (n, k) == (2, 2):
                assert literal_terminal == (0, 4, 1)
                assert shifted_terminal == (1, 4, 0)
                assert pword == (3, 2, 0)
            print('EXACT n=', n, 'k=', k, 'word=', pword,
                  'literal_terminal=', literal_terminal,
                  'zero_based_terminal=', shifted_terminal,
                  'all_three_internal_labelings_match=True', flush=True)

    print('SCOPE: exact finite comparison 0<=n<=9, 1<=k<=4; both word forms.')
    print('CONCLUSION: internal variant passes; source terminal convention differs.')
    print('No settlement of the published Conjecture 5 is claimed.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
