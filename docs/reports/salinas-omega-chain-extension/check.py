#!/usr/bin/env python3
"""Exact finite checks for Salinas Schmeis, arXiv:2610.29056v1, Q2.32.

Only integers, tuples and Boolean relations are used. Infinite orders are
given by formulas, but assertions are checked only on the reported windows.
This is experimental evidence, not a proof of the unbounded statement.
"""

from itertools import combinations, product


def order_axioms(points, le):
    table = [[le(x, y) for y in points] for x in points]
    for i in range(len(points)):
        assert table[i][i]
        for j in range(len(points)):
            assert not (i != j and table[i][j] and table[j][i])
            if table[i][j]:
                for k in range(len(points)):
                    assert not table[j][k] or table[i][k]


def induced_embedding(source, source_le, images, target_le):
    return len(set(images)) == len(source) and all(
        source_le(x, y) == target_le(images[i], images[j])
        for i, x in enumerate(source) for j, y in enumerate(source)
    )


def chain_prefix(sequence, le):
    assert induced_embedding(list(range(len(sequence))), lambda i, j: i <= j,
                             sequence, le)


def check_partition(points, parts, le):
    flattened = [x for c in parts for x in c]
    assert len(flattened) == len(set(flattened))
    assert set(flattened) == set(points)
    for c in parts:
        assert c
        chain_prefix(c, le)


def naturally_labeled_posets(n):
    """Every isomorphism type occurs: choose a linear extension as labeling."""
    pairs = list(combinations(range(n), 2))
    for bits in range(1 << len(pairs)):
        table = [[i == j for j in range(n)] for i in range(n)]
        for b, (i, j) in enumerate(pairs):
            table[i][j] = bool(bits & (1 << b))
        if all(not (table[i][j] and table[j][k]) or table[i][k]
               for i in range(n) for j in range(n) for k in range(n)):
            yield table


def merged_le(table, x, y):
    if x[0] == 'p':
        return y[0] == 't' or table[x[1]][y[1]]
    return y[0] == 't' and x[1] <= y[1]


def separated_le(table, x, y):
    if x[0] == 'p':
        return table[x[1]][y[1]]
    return y[0] == 'r' and x[1] == y[1] and x[2] <= y[2]


def has_copy(table, target, le):
    """Backtrack all induced injections, pruning both order directions."""
    n = len(table)
    images = []

    def extend(i):
        if i == n:
            return True
        for z in target:
            if z in images:
                continue
            if all(table[i][j] == le(z, y) and table[j][i] == le(y, z)
                   for j, y in enumerate(images)):
                images.append(z)
                if extend(i + 1):
                    return True
                images.pop()
        return False

    return extend(0)


def check_finite_posets():
    counts, greatest_counts = [], []
    assignments = 0
    obstruction_components = 0
    directed_tests = 0
    for n in range(1, 6):
        count = greatest_count = 0
        for table in naturally_labeled_posets(n):
            count += 1
            base = [('p', i) for i in range(n)]
            ble = lambda x, y: merged_le(table, x, y)
            # Theorem 2.5 fidelity: P + omega is countable and has no maxima.
            # Its partition is p_i followed by t_i,t_(n+i),t_(2n+i),... .
            tail = [('t', k) for k in range(6 * n)]
            bpoints = base + tail
            parts = [[('p', i)] + [('t', n * k + i) for k in range(6)]
                     for i in range(n)]
            check_partition(bpoints, parts, ble)
            order_axioms(bpoints, ble)
            assert induced_embedding(list(range(n)), lambda i, j: table[i][j],
                                     base, ble)
            for x in bpoints:
                successor = ('t', 0 if x[0] == 'p' else x[1] + 1)
                assert x != successor and ble(x, successor)
            # Explicit tail-escape witness used in the directedness argument.
            for i in range(n):
                for start in range(5):
                    for y in bpoints:
                        k = start + (y[1] + 1 if y[0] == 't' else 0)
                        z = ('t', n * k + i)
                        assert ble(parts[i][0], z) and ble(y, z)
                        directed_tests += 1
            maxima = [i for i in range(n)
                      if all(not table[i][j] or i == j for j in range(n))]
            tops = [i for i in range(n) if all(table[j][i] for j in range(n))]
            assert bool(tops) == (len(maxima) == 1)
            greatest_count += bool(tops)
            ale = lambda x, y: separated_le(table, x, y)
            apoints = base + [('r', m, k) for m in maxima for k in range(6)]
            order_axioms(apoints, ale)
            assert induced_embedding(list(range(n)), lambda i, j: table[i][j],
                                     base, ale)
            # All choices p <= assigned maximal m, rather than one convenient choice.
            choices = [[m for m in maxima if table[i][m]] for i in range(n)]
            for assignment in product(*choices):
                assignments += 1
                aparts = []
                for m in maxima:
                    fiber = [i for i in range(n) if assignment[i] == m]
                    assert fiber  # m itself is forced into this fiber.
                    width = len(fiber)
                    for residue, i in enumerate(fiber):
                        aparts.append([('p', i)] +
                                      [('r', m, k) for k in range(6)
                                       if k % width == residue])
                check_partition(apoints, aparts, ale)
            if not tops:
                for m in maxima:
                    component = [('p', i) for i in range(n) if table[i][m]]
                    component += [('r', m, k) for k in range(n)]
                    # n ray points suffice for an n-point copy after rank compression.
                    assert not has_copy(table, component, ale)
                    obstruction_components += 1
            else:
                t = tops[0]
                # Positive classification's map into the separated admissible extension.
                images = base + [('r', t, k) for k in range(6 * n)]
                assert induced_embedding(bpoints, ble, images, ale)
        counts.append(count)
        greatest_counts.append(greatest_count)
    assert counts == [1, 2, 7, 40, 357]
    print('FIDELITY Theorem 2.5: P + omega partitions for all naturally labeled posets')
    print('  n=1..5 poset counts:', counts, 'total:', sum(counts))
    print('  tail indices 0..6*n-1; strict successors also tested beyond the window')
    print('CLASSIFICATION greatest-element counts:', greatest_counts)
    print('  maximal-assignment partitions:', assignments)
    print('  no-copy lower-cone-plus-ray checks:', obstruction_components)
    print('  merged-tail escape witnesses:', directed_tests)


def check_connected_star_family():
    incompatible_pairs = 0
    forced_upper_triples = 0
    for m in range(2, 7):
        table = [[i == j or i == 0 for j in range(m + 1)]
                 for i in range(m + 1)]
        base = [('p', i) for i in range(m + 1)]
        ble = lambda x, y: merged_le(table, x, y)
        root = ('root',)

        def ale(x, y):
            return x == root or (y != root and x[1] == y[1] and x[2] <= y[2])

        for height in range(1, 9):
            aa = [root] + [('r', i, k) for i in range(m) for k in range(height)]
            bb = base + [('t', k) for k in range((m + 1) * height)]
            order_axioms(aa, ale)
            order_axioms(bb, ble)
            aimages = [root] + [('r', i, 0) for i in range(m)]
            assert induced_embedding(list(range(m + 1)), lambda i, j: table[i][j],
                                     aimages, ale)
            apart = [[('r', i, k) for k in range(height)] for i in range(m)]
            apart[0].insert(0, root)
            check_partition(aa, apart, ale)
            bpart = [[('p', 0), ('p', 1)] + [('t', m * k) for k in range(height)]]
            bpart += [[('p', i + 1)] + [('t', m * k + i) for k in range(height)]
                      for i in range(1, m)]
            bwindow = base + [('t', k) for k in range(m * height)]
            check_partition(bwindow, bpart, ble)
            for u, v in combinations(aa, 2):
                if not ale(u, v) and not ale(v, u):
                    assert u[0] == v[0] == 'r' and u[1] != v[1]
                    assert not any(ale(u, z) and ale(v, z) for z in aa)
                    incompatible_pairs += 1
            for u, v in product(bb, repeat=2):
                if not ble(u, v) and not ble(v, u):
                    assert u[0] == v[0] == 'p' and u[1] != 0 and v[1] != 0
                    for z in bb:
                        if u != z and ble(u, z):
                            assert z[0] == 't' and ble(v, z)
                            forced_upper_triples += 1
        print('COUNTEREXAMPLE connected star maxima=%d: ray heights 1..8 OK' % m)
    print('  separated-ray incomparable pairs:', incompatible_pairs)
    print('  merged-tail forced-upper-bound triples:', forced_upper_triples)


def main():
    # Empty P is not excluded by the source; omega is an admissible extension.
    chain_prefix([('t', k) for k in range(32)], lambda x, y: x[1] <= y[1])
    check_finite_posets()
    check_connected_star_family()
    print('SCOPE finite exact checks; infinite conclusion requires the preregistered proof')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
