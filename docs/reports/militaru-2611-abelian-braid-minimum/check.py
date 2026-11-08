#!/usr/bin/env python3
"""Exact finite checks for Militaru's Remark 3.4, finite abelian subfamily.

Standard library only. Enumerate ALL carrier permutations, identify the full
centralizer by an explicit group isomorphism, and check the orbit-kernel lower
bound over ALL subgroups and ALL possible kernel intersections in each case.
These finite checks support a general proof; they do not prove that theorem.
"""
from collections import Counter, deque
from itertools import permutations, product
from math import factorial, lcm, prod


def compose(f, g):
    return tuple(f[g[x]] for x in range(len(f)))


def perm_order(f):
    seen, answer = set(), 1
    for x in range(len(f)):
        if x in seen:
            continue
        y, length = x, 0
        while y not in seen:
            seen.add(y)
            length += 1
            y = f[y]
        answer = lcm(answer, length)
    return answer


def verify_braid(table):
    n = len(table)
    def b12(t):
        u, v = table[t[0]][t[1]]
        return u, v, t[2]
    def b23(t):
        u, v = table[t[1]][t[2]]
        return t[0], u, v
    for t in product(range(n), repeat=3):
        assert b12(b23(b12(t))) == b23(b12(b23(t))), t
    assert len({table[x][y] for x in range(n) for y in range(n)}) == n * n
    for x in range(n):
        assert len({table[x][y][0] for y in range(n)}) == n
    for y in range(n):
        assert len({table[x][y][1] for x in range(n)}) == n


def full_automorphisms(table):
    n = len(table)
    answer, scanned = set(), 0
    for f in permutations(range(n)):
        scanned += 1
        if all((f[table[x][y][0]], f[table[x][y][1]]) == table[f[x]][f[y]]
               for x in range(n) for y in range(n)):
            answer.add(f)
    assert scanned == factorial(n)
    return answer, scanned


def fidelity():
    # Example 3.3 / Remark 3.4: R(x,y)=(y,x+1) on Z/nZ,
    # and its FULL automorphisms are precisely the cyclic translations.
    for n in range(2, 9):
        table = [[(y, (x + 1) % n) for y in range(n)] for x in range(n)]
        verify_braid(table)
        aut, scanned = full_automorphisms(table)
        translations = {tuple((x + c) % n for x in range(n)) for c in range(n)}
        assert aut == translations
        gen = tuple((x + 1) % n for x in range(n))
        assert perm_order(gen) == n
        print(f'FIDELITY cyclic n={n}: permutations={scanned}; braid, bijective, '
              f'both nondegenerate; full Aut=C{n}')


def prime_divisors(n):
    result, p = [], 2
    while p * p <= n:
        if n % p == 0:
            result.append(p)
            while n % p == 0:
                n //= p
        p += 1
    if n > 1:
        result.append(n)
    return result


def exact_power_log(n, p):
    answer = 0
    while n > 1:
        assert n % p == 0
        n //= p
        answer += 1
    return answer


def weight(elements, add, zero):
    """Compute sum of primary orders from the intrinsic p-torsion filtration."""
    def scale(k, x):
        y = zero
        for _ in range(k):
            y = add(y, x)
        return y
    total = 0
    for p in prime_divisors(len(elements)):
        exponent, bound = 0, 1
        while bound <= len(elements):
            exponent += 1
            bound *= p
        hp = {x for x in elements if scale(bound, x) == zero}
        last = hp
        k = 1
        while len(last) > 1:
            torsion = {x for x in last if scale(p, x) == zero}
            d = exact_power_log(len(torsion), p)
            total += (p if k == 1 else p ** k - p ** (k - 1)) * d
            last = {scale(p, x) for x in last}
            k += 1
    return total


def lower_bound(qs):
    elements = tuple(product(*(range(q) for q in qs)))
    zero = (0,) * len(qs)
    def add(x, y):
        return tuple((a + b) % q for a, b, q in zip(x, y, qs))
    def extend(h, x):
        result, frontier = set(h), list(h)
        while frontier:
            y = add(frontier.pop(), x)
            if y not in result:
                result.add(y)
                frontier.append(y)
        return frozenset(result)
    trivial = frozenset([zero])
    subgroups, frontier = {trivial}, [trivial]
    while frontier:
        h = frontier.pop()
        for x in elements:
            k = extend(h, x)
            if k not in subgroups:
                subgroups.add(k)
                frontier.append(k)
    whole = frozenset(elements)
    weights = {h: weight(h, add, zero) for h in subgroups}
    assert weights[whole] == sum(qs)
    quotient_weight, index = {}, {}
    for h in subgroups:
        representative = {x: min(add(x, y) for y in h) for x in elements}
        quotient = frozenset(representative.values())
        quotient_add = lambda x, y: representative[add(x, y)]
        quotient_weight[h] = weight(quotient, quotient_add, representative[zero])
        index[h] = len(quotient)
        assert index[h] * len(h) == len(elements)
        assert quotient_weight[h] <= index[h]
    for h in subgroups:
        for k in subgroups:
            if h <= k:
                assert weights[h] <= weights[k]
            # A/(H cap K) embeds in A/H x A/K.
            intersection = h & k
            assert quotient_weight[intersection] <= quotient_weight[h] + quotient_weight[k]
            # Independently compute the product weight, checking additivity.
            cartesian = tuple(product(h, k))
            cart_add = lambda x, y: (add(x[0], y[0]), add(x[1], y[1]))
            assert weight(cartesian, cart_add, (zero, zero)) == weights[h] + weights[k]
    # Dynamic programming over kernel intersections checks arbitrary families
    # of coset orbits, without imposing a bound on the number of orbits.
    costs, queue = {whole: 0}, deque([whole])
    while queue:
        kernel = queue.popleft()
        for h in subgroups:
            new_kernel, cost = kernel & h, costs[kernel] + index[h]
            if new_kernel not in costs or cost < costs[new_kernel]:
                costs[new_kernel] = cost
                queue.append(new_kernel)
    for kernel, cost in costs.items():
        assert quotient_weight[kernel] <= cost
    assert costs[trivial] == sum(qs)
    return len(subgroups), costs[trivial]


def triangular_case(label, qs):
    assert all(len(prime_divisors(q)) == 1 for q in qs)
    points = tuple((j, b) for j, q in enumerate(qs) for b in range(q))
    address = {x: i for i, x in enumerate(points)}
    n = len(points)
    left = tuple(tuple(address[j, (b + (i <= j)) % qs[j]] for j, b in points)
                 for i in range(len(qs)))
    for i, lam in enumerate(left):
        for j in range(len(qs)):
            assert any(lam[address[j, b]] != address[j, b] for b in range(qs[j])) == (i <= j)
    assert len(set(left)) == len(qs)
    table = [[(left[points[x][0]][y], x) for y in range(n)] for x in range(n)]
    verify_braid(table)
    aut, scanned = full_automorphisms(table)
    # Explicit isomorphism c -> (j,b) |-> (j,b+c_j), with a tested inverse
    # and the multiplication table; order distributions alone do not suffice.
    coords = tuple(product(*(range(q) for q in qs)))
    phi = {c: tuple(address[j, (b + c[j]) % qs[j]] for j, b in points) for c in coords}
    assert len(set(phi.values())) == len(coords)
    assert set(phi.values()) == aut
    for c, f in phi.items():
        recovered = tuple(points[f[address[j, 0]]][1] for j in range(len(qs)))
        assert recovered == c
        for d, g in phi.items():
            c_plus_d = tuple((a + b) % q for a, b, q in zip(c, d, qs))
            assert compose(f, g) == phi[c_plus_d]
    distribution = dict(sorted(Counter(map(perm_order, aut)).items()))
    subgroups, minimum = lower_bound(qs)
    if label == 'C6':
        # CRT identifies the stated C6 example with its PRIMARY factors C2 x C3.
        crt = {a: (a % 2, a % 3) for a in range(6)}
        assert set(crt.values()) == set(coords)
        for a, b in product(range(6), repeat=2):
            assert compose(phi[crt[a]], phi[crt[b]]) == phi[crt[(a + b) % 6]]
    print(f'CASE {label}: primary={qs}; |X|={n}; permutations={scanned}; '
          f'|Aut|={len(aut)}; explicit product isomorphism OK; orders={distribution}; '
          f'all {subgroups} subgroups checked; faithful degree={minimum}')


def main():
    fidelity()
    for label, qs in [('C2^2', (2, 2)), ('C2 x C4', (2, 4)),
                      ('C4 x C2', (4, 2)), ('C3^2', (3, 3)), ('C6', (2, 3)),
                      ('C2^3', (2, 2, 2)), ('C4^2', (4, 4)),
                      ('C2^2 x C3', (2, 2, 3))]:
        triangular_case(label, qs)
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
