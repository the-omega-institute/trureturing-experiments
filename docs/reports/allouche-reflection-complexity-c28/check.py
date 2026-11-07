#!/usr/bin/env python3
"""Exact finite checks for ACLSS Conjecture 28; no universal settlement is asserted."""
from itertools import product


def windows(word, n):
    return {tuple(word[i:i + n]) for i in range(len(word) - n + 1)}


def orbit(word):
    return min(word, word[::-1])


def morphism(word):
    return tuple(b for a in word for b in ((0, 1) if a == 0 else (1, 0)))


def tm_factors(n):
    # Every length-n factor meets at most two aligned mu^k blocks, 2^k >= n.
    # All four binary pairs occur in Thue--Morse; their mu^k-images occur too.
    images = [p for p in product(range(2), repeat=2)]
    block_size = 1
    while block_size < n:
        images = [morphism(p) for p in images]
        block_size *= 2
    return set().union(*(windows(w, n) for w in images))


def fidelity():
    expected = [1, 2, 3, 4, 6, 6, 10, 10, 13, 12, 16, 16, 20, 20, 22]
    actual = [len({orbit(w) for w in tm_factors(n)}) for n in range(len(expected))]
    assert actual == expected, (actual, expected)
    assert tm_factors(2) == set(product(range(2), repeat=2))
    print('Fidelity: source Example 4 / equation (rt), Thue--Morse r(0..14) =', actual)


def check_word(prefix, period, max_n):
    m, t = len(prefix), len(period)

    def x(i):
        return prefix[i] if i < m else period[(i - m) % t]

    def z(i):
        return -1 if i == 0 else x(abs(i) - 1)

    def language(n):
        # Starts >= m repeat modulo t; these are ALL factors of prefix.period^omega.
        return {tuple(x(i + j) for j in range(n)) for i in range(m + t)}

    def mirror_language(n):
        # Includes all right-tail starts, left-tail starts, and all marker positions.
        return {tuple(z(i + j) for j in range(n))
                for i in range(-n - m - t + 1, m + t + 1)}

    data = []
    for n in range(max_n + 2):
        lang = language(n)
        sym = lang | {w[::-1] for w in lang}
        classes = {orbit(w) for w in lang}
        pal = sum(w == w[::-1] for w in sym)
        mir = mirror_language(n)
        marked = {w for w in mir if -1 in w}
        assert {w for w in mir if -1 not in w} == sym
        assert len(marked) == n
        assert {w.index(-1) for w in marked} == set(range(n))
        assert len(mir) == len(sym) + n
        assert 2 * len(classes) == len(sym) + pal
        data.append((sym, classes, pal, mir))

    negative_steps = 0
    non_reversal_closed = False
    for n in range(max_n + 1):
        sym, vertices, pn, mir = data[n]
        sym_next, edges, pn_next, mir_next = data[n + 1]
        parent = {v: v for v in vertices}

        def root(v):
            while parent[v] != v:
                v = parent[v]
            return v

        for w in edges:
            a, b = orbit(w[:-1]), orbit(w[1:])
            assert a in vertices and b in vertices
            assert {a, b} == {orbit(w[::-1][:-1]), orbit(w[::-1][1:])}
            if w == w[::-1]:
                assert a == b
            else:
                parent[root(a)] = root(b)
        # Connected even after deleting ALL palindromic edges, retaining other loops.
        assert len({root(v) for v in vertices}) == 1
        assert len(edges) - pn_next >= len(vertices) - 1
        du = len(sym_next) - len(sym)
        dq = len(edges) - len(vertices)
        assert 2 * dq == du + pn_next - pn
        assert pn + pn_next <= du + 2
        assert -1 <= dq <= du + 1
        assert du == len(mir_next) - len(mir) - 1
        negative_steps += dq < 0
        non_reversal_closed |= language(n) != sym
    return negative_steps, non_reversal_closed


def main():
    fidelity()
    max_n = 16
    cases = comparisons = negative_steps = nonclosed = 0
    # Complete Cartesian products in these declared ranges; no random sampling.
    for alphabet_size, prefix_max, period_max in [(1, 3, 5), (2, 3, 5), (3, 2, 3)]:
        alphabet = range(alphabet_size)
        group = 0
        for prefix_len in range(prefix_max + 1):
            for prefix in product(alphabet, repeat=prefix_len):
                for period_len in range(1, period_max + 1):
                    for period in product(alphabet, repeat=period_len):
                        neg, nc = check_word(prefix, period, max_n)
                        negative_steps += neg
                        nonclosed += nc
                        cases += 1
                        group += 1
                        comparisons += max_n + 1
        print(f'Exhaustive family: alphabet={alphabet_size}, prefix length=0..{prefix_max}, '
              f'period length=1..{period_max}: {group} presentations')
    assert cases == 1457
    assert nonclosed > 0
    assert negative_steps > 0
    print(f'Checked {cases} ultimately periodic presentations, {comparisons} consecutive-length '
          f'comparisons (n=0..{max_n}); {nonclosed} presentations lack reversal closure in the range.')
    print(f'Observed {negative_steps} negative reflection-complexity steps; integer subtraction preserved.')
    print('Verified: marked-mirror factor bijection/count; orbit identity; quotient endpoints; '
          'palindromic loops; connectivity after loop deletion; palindromic sum and difference inequalities.')
    print('Scope: finite exact checks; neither Cassaigne\'s infinite-language theorem nor the universal '
          'bounded-difference conjecture is proved by this program.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
