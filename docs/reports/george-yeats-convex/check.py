#!/usr/bin/env python3
"""Exact checks for the George–Yeats s=1 specialization; standard library only.

Polynomials are Counters from exponent tuples to integer coefficients. Posets
are tuples of bit masks: row a is the reflexive upper set of a. No floating
point arithmetic, random sampling, external solvers or source proof is used.
"""

from collections import Counter
from functools import lru_cache
from itertools import combinations
from math import comb


def bits(mask):
    while mask:
        low = mask & -mask
        yield low.bit_length() - 1
        mask ^= low


def induced(p, vertices):
    return tuple(sum(1 << j for j, b in enumerate(vertices) if p[a] >> b & 1)
                 for a in vertices)


def combine(a, b, series):
    n, m = len(a), len(b)
    upper = ((1 << m) - 1) << n if series else 0
    return tuple(row | upper for row in a) + tuple(row << n for row in b)


def normalize(op, a, b):
    children = (a[1] if a[0] == op else (a,)) + (b[1] if b[0] == op else (b,))
    return (op, tuple(sorted(children)) if op == 'p' else children)


@lru_cache(None)
def realize(tree):
    if tree == ('a',):
        return (1,)
    p = ()
    for child in tree[1]:
        p = combine(p, realize(child), tree[0] == 's')
    return p


def components(p, comparable):
    n = len(p)
    unseen, answer = set(range(n)), []
    while unseen:
        todo = [min(unseen)]
        unseen.remove(todo[0])
        comp = []
        while todo:
            a = todo.pop()
            comp.append(a)
            adjacent = [b for b in unseen
                        if bool((p[a] >> b & 1) or (p[b] >> a & 1)) == comparable]
            unseen.difference_update(adjacent)
            todo.extend(adjacent)
        answer.append(tuple(sorted(comp)))
    return answer


def decompose(p):
    """An isomorphism invariant obtained directly from the order relation."""
    if len(p) == 1:
        return ('a',)
    cs = components(p, True)
    if len(cs) > 1:
        sub = [decompose(induced(p, c)) for c in cs]
        return None if None in sub else ('p', tuple(sorted(sub)))
    cs = components(p, False)
    if len(cs) < 2:
        return None
    all_cs = tuple(cs)
    cs.sort(key=lambda c: sum(bool(p[b[0]] >> c[0] & 1) for b in all_cs))
    for i, a in enumerate(cs):
        for b in cs[i + 1:]:
            assert all(p[x] >> y & 1 for x in a for y in b)
    sub = [decompose(induced(p, c)) for c in cs]
    return None if None in sub else ('s', tuple(sub))


def n_free(p):
    # The comparability graph of an induced N is the four-vertex path.
    for subset in combinations(range(len(p)), 4):
        degrees = [sum(bool((p[a] >> b & 1) or (p[b] >> a & 1))
                       for b in subset if b != a) for a in subset]
        if sorted(degrees) == [1, 1, 2, 2]:
            return False
    return True


def statistics(p):
    n = len(p)
    downs = tuple(sum(1 << a for a in range(n) if p[a] >> b & 1) for b in range(n))
    rows = []
    for x in range(1 << n):
        d = u = mins = maxs = 0
        for a in bits(x):
            d |= downs[a]
            u |= p[a]
            if downs[a] & x == 1 << a:
                mins |= 1 << a
            if p[a] & x == 1 << a:
                maxs |= 1 << a
        if d & u != x:
            continue
        # Independent literal interval criterion for source fidelity.
        assert all((p[a] & downs[b]) & ~x == 0 for a in bits(x) for b in bits(x))
        literal_d = literal_u = 0
        for a in bits(maxs):
            literal_d |= downs[a]
        for a in bits(mins):
            literal_u |= p[a]
        assert (literal_d, literal_u) == (d, u)
        rows.append((x, mins.bit_count(), maxs.bit_count(), d, u))
    return rows


def polynomials(p):
    phi, psi = Counter(), Counter()
    data = statistics(p)
    for x, s, t, d, u in data:
        exponent = (s, t, d.bit_count(), u.bit_count(), x.bit_count())
        phi[exponent] += 1
        psi[exponent[1:]] += 1
    return phi, psi, data


def multiply(a, b):
    out = Counter()
    for u, c in a.items():
        for v, d in b.items():
            out[tuple(x + y for x, y in zip(u, v))] += c * d
    return out


def fidelity():
    # Definition 2.1 / Example 2.2 and the full Example 2.3 diamond formula.
    assert polynomials((1,))[0] == Counter({(0, 0, 0, 0, 0): 1, (1, 1, 1, 1, 1): 1})
    diamond = (15, 10, 12, 8)
    bracket = {
        (0, 0, 3, 3, 3): 1, (1, 0, 3, 2, 2): 1,
        (0, 1, 2, 3, 2): 1, (1, 1, 2, 2, 1): 1,
        (0, 0, 3, 1, 1): 2, (0, 0, 1, 3, 1): 2,
        (0, 0, 3, 0, 0): 1, (0, 0, 0, 3, 0): 1,
        (0, 0, 1, 1, 0): 2,
    }
    expected = Counter({(0, 0, 0, 0, 0): 1})
    for e, c in bracket.items():
        expected[tuple(a + 1 for a in e)] += c
    actual = polynomials(diamond)[0]
    assert actual == expected
    for n in range(1, 9):
        antichain = tuple(1 << i for i in range(n))
        assert polynomials(antichain)[0] == Counter({(k,) * 5: comb(n, k)
                                                   for k in range(n + 1)})
    print('FIDELITY: singleton, full five-variable diamond formula (13 convex sets),')
    print('          antichain formula (1+stxyz)^n for 1 <= n <= 8: OK')


def generate(bound):
    levels = {1: {('a',)}}
    for n in range(2, bound + 1):
        levels[n] = set()
        for k in range(1, n):
            for a in levels[k]:
                for b in levels[n - k]:
                    for op in ('p', 's'):
                        levels[n].add(normalize(op, a, b))
    return levels


def independent_census(levels):
    # Every finite poset has a linear extension; all natural labellings suffice.
    total = 0
    for n in range(1, 7):
        pairs = list(combinations(range(n), 2))
        seen = set()
        valid = 0
        for mask in range(1 << len(pairs)):
            p = [1 << i for i in range(n)]
            for bit, (a, b) in enumerate(pairs):
                if mask >> bit & 1:
                    p[a] |= 1 << b
            if any(p[b] & ~p[a] for a in range(n) for b in bits(p[a])):
                continue
            valid += 1
            signature = decompose(p)
            assert (signature is not None) == n_free(p)
            if signature is not None:
                seen.add(signature)
        assert seen == levels[n]
        total += valid
        print(f'INDEPENDENT_CENSUS: n={n}, naturally-labelled posets={valid}, SP classes={len(seen)}: OK')
    return total


def check_structure(p, psi, data):
    n = len(p)
    assert psi[(0, 0, 0, 0)] == 1
    assert all(all(a > 0 for a in e) for e in psi if any(e))
    assert max(e[2] for e in psi) == max(e[3] for e in psi) == n
    if len(components(p, True)) == 1:
        assert any(e[0] == 1 and e[2] == n and e[1] == e[3] for e in psi)

    down_poly = Counter()
    for (t, d, u, z), c in psi.items():
        if d == z:
            down_poly[d] += c
    ideals = [(x, x.bit_count()) for x, _, _, d, _ in data if x == d]
    assert down_poly == Counter(k for _, k in ideals)
    cuts = [k for k in range(n + 1) if down_poly[k] == 1]
    cut_sets = []
    for k in range(n + 1):
        at_k = [x for x, size in ideals if size == k]
        ordinal = [x for x in at_k if all(p[a] >> b & 1
                   for a in bits(x) for b in range(n) if not x >> b & 1)]
        assert (len(at_k) == 1) == bool(ordinal)
        if len(at_k) == 1:
            cut_sets.append(at_k[0])
    assert all(a & ~b == 0 for a, b in zip(cut_sets, cut_sets[1:]))
    assert cuts[0] == 0 and cuts[-1] == n
    for i in range(len(cuts) - 1):
        lo, hi = cuts[i:i + 2]
        recovered = Counter({(0, 0, 0, 0): 1})
        for (t, d, u, z), c in psi.items():
            if lo < d <= hi and n - hi < u <= n - lo:
                recovered[t, d - lo, u - n + hi, z] += c
        block = tuple(bits(cut_sets[i + 1] & ~cut_sets[i]))
        assert recovered == polynomials(induced(p, block))[1]


def main():
    fidelity()
    levels = generate(8)
    independent_census(levels)
    saved, checked = {}, 0
    for n in range(1, 9):
        by_polynomial = {}
        connected = 0
        for tree in sorted(levels[n]):
            p = realize(tree)
            assert decompose(p) == tree
            _, psi, data = polynomials(p)
            check_structure(p, psi, data)
            key = tuple(sorted(psi.items()))
            if key in by_polynomial:
                raise AssertionError(('SPECIALIZATION_COLLISION', tree, by_polynomial[key]))
            by_polynomial[key] = tree
            saved[tree] = psi
            connected += len(components(p, True)) == 1
            checked += 1
        print(f'SPECIALIZATION: n={n}, SP isomorphism classes={len(levels[n])}, connected={connected}, collisions=0: OK')
    products = 0
    for n in range(2, 8):
        for k in range(1, n):
            for a in levels[k]:
                for b in levels[n - k]:
                    assert saved[normalize('p', a, b)] == multiply(saved[a], saved[b])
                    products += 1
    print(f'STRUCTURE: {checked} classes; degree/support, connected t^1 y^n witness,')
    print('           downset filter, unique-downset cuts and two-coordinate block extraction: OK')
    print(f'DISJOINT_UNION: {products} ordered summand pairs with total size <= 7: OK')
    print('SCOPE: exact finite verification through 8 elements; no unbounded or Lean proof is asserted.')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
