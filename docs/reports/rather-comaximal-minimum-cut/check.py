#!/usr/bin/env python3
"""Exact residue-graph fidelity and finite minimum-cut checks for Rather 8.1."""
import itertools
import json
import math


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def factor(n):
    primes = []
    squarefree = True
    d = 2
    while d * d <= n:
        exponent = 0
        while n % d == 0:
            n //= d
            exponent += 1
        if exponent:
            primes.append(d)
            squarefree = squarefree and exponent == 1
        d += 1
    if n > 1:
        primes.append(n)
    return primes, squarefree


def graph(n):
    vertices = [x for x in range(1, n) if math.gcd(x, n) > 1]
    neighbors = [0] * len(vertices)
    for i, x in enumerate(vertices):
        for j in range(i):
            if math.gcd(math.gcd(x, vertices[j]), n) == 1:
                neighbors[i] |= 1 << j
                neighbors[j] |= 1 << i
    return vertices, neighbors


def preconnected(neighbors, remaining):
    if remaining == 0:
        return True
    reached = remaining & -remaining
    frontier = reached
    while frontier:
        bit = frontier & -frontier
        frontier ^= bit
        new = neighbors[bit.bit_length() - 1] & remaining & ~reached
        reached |= new
        frontier |= new
    return reached == remaining


def support(x, primes):
    return sum(1 << i for i, p in enumerate(primes) if x % p == 0)


def fidelity(n, expected_degree):
    primes, squarefree = factor(n)
    require(squarefree and len(primes) >= 2, 'source modulus domain')
    vertices, neighbors = graph(n)
    minimum_degree = min(row.bit_count() for row in neighbors)
    require(minimum_degree == expected_degree, f'source minimum-degree table n={n}')
    counts = {}
    for x in vertices:
        s = support(x, primes)
        counts[s] = counts.get(s, 0) + 1
    for s in range(1, (1 << len(primes)) - 1):
        predicted = math.prod(p - 1 for i, p in enumerate(primes) if not s & (1 << i))
        require(counts[s] == predicted, f'Proposition 3.2 n={n}, S={s}')
    # Theorem 3.3, checked against the gcd-built graph, including its diagonal.
    supports = [support(x, primes) for x in vertices]
    for i in range(len(vertices)):
        for j in range(len(vertices)):
            predicted = i != j and supports[i] & supports[j] == 0
            require(bool(neighbors[i] & (1 << j)) == predicted, f'adjacency n={n}')
    distinguished = [i for i, s in enumerate(supports) if s == 1 << (len(primes) - 1)]
    require(len(distinguished) == expected_degree, f'Lemma 5.1 size n={n}')
    deleted = sum(1 << i for i in distinguished)
    remaining = ((1 << len(vertices)) - 1) ^ deleted
    require(not preconnected(neighbors, remaining), f'Lemma 5.1 disconnection n={n}')
    return {'n': n, 'vertices': len(vertices), 'support_classes': len(counts),
            'minimum_degree': minimum_degree,
            'distinguished_cut': [vertices[i] for i in distinguished]}


def enumerate_cuts(n):
    primes, _ = factor(n)
    vertices, neighbors = graph(n)
    delta = math.prod(p - 1 for p in primes[:-1])
    expected = tuple(x for x in vertices if support(x, primes) == 1 << (len(primes) - 1))
    all_vertices = (1 << len(vertices)) - 1
    checks = []
    cuts = []
    for k in range(delta + 1):
        checked = 0
        for c in itertools.combinations(range(len(vertices)), k):
            removed = sum(1 << i for i in c)
            if not preconnected(neighbors, all_vertices ^ removed):
                cuts.append(tuple(vertices[i] for i in c))
                require(k == delta, f'cut smaller than delta n={n}: {c}')
            checked += 1
        require(checked == math.comb(len(vertices), k), 'enumeration count')
        checks.append(checked)
    require(cuts == [expected], f'unique minimum-cut classification n={n}: {cuts}')
    return {'n': n, 'primes': primes, 'vertices': len(vertices), 'delta': delta,
            'checked_deletions_by_size': checks, 'minimum_cuts': cuts}


def main():
    # Table 2 in the TeX: exact kappa and minimum degree are the listed integers.
    source = [fidelity(n, d) for n, d in
              [(6, 1), (30, 2), (42, 2), (70, 4), (210, 8), (2310, 48)]]
    cases = []
    for n in range(2, 71):
        primes, squarefree = factor(n)
        if squarefree and len(primes) >= 2:
            cases.append(enumerate_cuts(n))
    seed_cases = {6, 10, 15, 30, 35, 42, 70}
    seed_delta_deletions = sum(x['checked_deletions_by_size'][-1]
                               for x in cases if x['n'] in seed_cases)
    require(seed_delta_deletions == 149844, 'search-seat count independently reproduced')
    result = {
        'source': 'arXiv:2605.01565v1, Proposition 3.2, Theorem 3.3, Lemma 5.1, Table 2',
        'arithmetic': 'exact Python integers and gcd; no floating point',
        'fidelity': source,
        'classification_scope': 'every squarefree composite integer 2 <= n <= 70; all deletions of size <= delta',
        'case_count': len(cases),
        'total_deletions_checked': sum(sum(x['checked_deletions_by_size']) for x in cases),
        'search_seat_seven_case_delta_deletions': seed_delta_deletions,
        'cases': cases,
        'boundary': 'finite checks are fidelity helpers and do not prove the universal classification',
    }
    print(json.dumps(result, indent=2, sort_keys=True))
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
