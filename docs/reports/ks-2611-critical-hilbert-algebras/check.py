#!/usr/bin/env python3
"""Exact exhaustive tests of the proposed finite Ho criticality criterion.

A finite algebra embeds into a product of proper subalgebras iff the
homomorphisms into proper subalgebras jointly separate its elements. Every
such map, followed by inclusion, is a nonsurjective endomorphism; conversely
an endomorphism's image is a subalgebra. The direct test below enumerates
these maps from the operation table, independently of the fold criterion.
Only Python integers, sets, tuples, and exact equality are used.
"""
from itertools import combinations, permutations, product
from math import factorial
import argparse


def subsets(s):
    t = s
    while True:
        yield t
        if not t:
            break
        t = (t - 1) & s


def posets(n):
    """One representative of each poset isomorphism class, without a library.

    Every finite poset has a linear extension. Enumerate all strict relations
    compatible with 0,...,n-1, keep precisely the transitive ones, then use all
    permutations to remove isomorphic duplicates.
    """
    pairs = list(combinations(range(n), 2))
    seen = set()
    for bits in range(1 << len(pairs)):
        up = [1 << i for i in range(n)]
        for k, (i, j) in enumerate(pairs):
            if bits >> k & 1:
                up[i] |= 1 << j
        if any(up[j] & ~up[i] for i in range(n) for j in range(n)
               if up[i] >> j & 1):
            continue
        encodings = []
        for perm in permutations(range(n)):
            encodings.append(tuple(sum(1 << perm[j] for j in range(n)
                                       if up[i] >> j & 1)
                                   for i in sorted(range(n), key=perm.__getitem__)))
        key = min(encodings, default=())
        if key not in seen:
            seen.add(key)
            yield tuple(up)


def maximal(up):
    return sum(1 << i for i, row in enumerate(up) if row == 1 << i)


def down(up, c):
    return sum(1 << i for i, row in enumerate(up) if row & c)


def antichain(up, c):
    return all((up[i] & c) == 1 << i for i in range(len(up)) if c >> i & 1)


def families(up):
    n = len(up)
    m = maximal(up)
    forced = set(subsets(m)) | {1 << i for i in range(n)}
    optional = sorted((c for c in range(1 << n)
                       if c not in forced and antichain(up, c)),
                      key=lambda c: (c.bit_count(), c))

    def visit(k, selected):
        if k == len(optional):
            yield tuple(sorted(selected))
            return
        c = optional[k]
        yield from visit(k + 1, selected)
        if all(c ^ (1 << i) in selected for i in range(n) if c >> i & 1):
            selected.add(c)
            yield from visit(k + 1, selected)
            selected.remove(c)
    yield from visit(0, forced)


def algebra(up, s):
    index = {c: i for i, c in enumerate(s)}
    table = tuple(tuple(index[b & ~down(up, c)] for b in s) for c in s)
    return table, index[0], index[maximal(up)]


def endomorphisms(table, one, zero):
    """Enumerate ALL constant- and implication-preserving maps by closure.

    Assign a source generator any possible target. Each assigned pair forces
    the image of its product. A conflict discards the branch. Branching on
    the next unassigned source point covers every homomorphism exactly once.
    """
    size = len(table)

    def extend(old, seeds):
        images = old.copy()
        queue = list(seeds)
        while queue:
            a, target = queue.pop()
            if images[a] >= 0:
                if images[a] != target:
                    return None
                continue
            images[a] = target
            for b, other in enumerate(images):
                if other >= 0:
                    queue.append((table[a][b], table[target][other]))
                    queue.append((table[b][a], table[other][target]))
        return images

    def visit(images):
        if -1 not in images:
            yield tuple(images)
            return
        a = images.index(-1)
        # Search near the identity first; this changes no accepted map.
        for target in (a, *(x for x in range(size) if x != a)):
            extended = extend(images, [(a, target)])
            if extended is not None:
                yield from visit(extended)
    images = extend([-1] * size, [(one, one), (zero, zero)])
    if images is not None:
        yield from visit(images)


def direct_critical(table, one, zero):
    """Search for a subdirect embedding into proper image subalgebras.

    The product of all proper images suffices: its kernel is the intersection
    of the coordinate kernels. Once a selected subfamily already separates
    every pair, its product is an embedding and enumeration can stop. A
    critical verdict always exhausts every endomorphism. Surjectivity onto
    each image is automatic. The empty product is permitted, so the
    one-element algebra is not critical.
    """
    size = len(table)
    inseparable = set(combinations(range(size), 2))
    proper = []
    for f in endomorphisms(table, one, zero):
        image = set(f)
        if len(image) < size:
            assert one in image and zero in image
            assert all(table[a][b] in image for a in image for b in image)
            assert all(f[table[a][b]] == table[f[a]][f[b]]
                       for a in range(size) for b in range(size))
            proper.append(f)
            inseparable = {(a, b) for a, b in inseparable if f[a] == f[b]}
            if not inseparable and size > 5:
                break
    if size <= 5:
        brute = []
        for f in product(range(size), repeat=size):
            if f[one] == one and f[zero] == zero and all(
                    f[table[a][b]] == table[f[a]][f[b]]
                    for a in range(size) for b in range(size)):
                brute.append(f)
        assert set(endomorphisms(table, one, zero)) == set(brute)
    return bool(inseparable), proper


def folds(up, s, d):
    m = maximal(up)
    ms = tuple(i for i in range(len(up)) if m >> i & 1)
    active = up[d] & m
    for values in product(ms, repeat=len(ms)):
        p = dict(zip(ms, values))
        if any(p[p[x]] != p[x] for x in ms):
            continue
        if any(p[x] != x for x in ms if active >> x & 1):
            continue
        if all(p[x] == x for x in ms):
            continue
        yield p


def transform(c, m, p):
    return (c & ~m) | sum(1 << x for x in p if c >> p[x] & 1)


def criterion(up, s):
    m = maximal(up)
    ds = [i for i in range(len(up)) if not (m >> i & 1)]
    if not ds:
        return m.bit_count() == 1
    least = [i for i in ds if all(up[i] >> j & 1 for j in ds)]
    if not least:
        return False
    faces = [c for c in s if not any(c != b and c & ~b == 0 for b in s)]
    for p in folds(up, s, least[0]):
        if all(transform(c, m, p) in s for c in faces):
            return False
    return True


def compose(f, g):
    return tuple(f[g[i]] for i in range(len(g)))


def power(f, exponent):
    result = tuple(range(len(f)))
    while exponent:
        if exponent & 1:
            result = compose(f, result)
        f = compose(f, f)
        exponent //= 2
    return result


def normalization_checks(up, s, table, one, proper):
    m = maximal(up)
    ds = [i for i in range(len(up)) if not (m >> i & 1)]
    least = [i for i in ds if all(up[i] >> j & 1 for j in ds)]
    if not least:
        return 0
    d = least[0]
    di = s.index(1 << d)
    index = {c: i for i, c in enumerate(s)}
    count = 0
    for f in proper:
        if f[di] == one:
            continue
        assert len({(f[i], c & m) for i, c in enumerate(s)}) == len(s)
        g = power(f, factorial(len(s)))
        assert compose(g, g) == g
        assert all(s[g[i]] == c for i, c in enumerate(s) if c & m == 0)
        ms = [i for i in range(len(up)) if m >> i & 1]
        p = {}
        for t in ms:
            matches = [a for a in ms if s[g[index[1 << a]]] >> t & 1]
            assert len(matches) == 1
            p[t] = matches[0]
        assert all(p[p[t]] == p[t] for t in ms)
        assert all(p[t] == t for t in ms if up[d] >> t & 1)
        assert all(s[g[i]] == transform(c, m, p) for i, c in enumerate(s))
        assert any(p[t] != t for t in ms)
        count += 1
    return count


def identities(up, s, table, one, zero):
    size = len(s)
    for x in range(size):
        assert table[x][x] == one and table[one][x] == x
        assert table[zero][x] == one
        dd = table[table[x][zero]][zero]
        assert s[dd] == s[x] & maximal(up)
        assert s[table[dd][x]] == s[x] & ~maximal(up)
        for y in range(size):
            assert table[table[x][y]][table[table[y][x]][x]] == \
                   table[table[y][x]][table[table[x][y]][y]]
            for z in range(size):
                assert table[x][table[y][z]] == table[table[x][y]][table[x][z]]


def source_example(n):
    # Indices m=0,a=1,b=2,c=3,e=4,s_i=5+i.
    count = n + 4
    up = [1 << i for i in range(count)]
    up[2] |= 1 << 1
    up[4] |= 1 << 3
    for i in range(n - 1):
        up[5 + i] |= (1 << 1) | (1 << 2) | (1 << 3) | (1 << 4)
        up[5 + i] |= sum(1 << (5 + j) for j in range(i))
    m = maximal(up)
    s = tuple(sorted(set(subsets(m)) | {1 << i for i in range(count)} |
                     {(1 << 1) | (1 << 4), (1 << 3) | (1 << 2)}))
    return tuple(up), s


def source_checks(n):
    up, s = source_example(n)
    table, one, zero = algebra(up, s)
    assert len(s) == n + 11  # Source's displayed cardinality.
    assert all(antichain(up, c) for c in s)
    assert all(b in s for c in s for b in subsets(c))
    identities(up, s, table, one, zero)
    def h(i):
        above = [j for j in range(len(up)) if j != i and up[i] >> j & 1]
        return 1 + max((h(j) for j in above), default=0)
    assert max(map(h, range(len(up)))) == n + 1
    # Source Proposition 4.1: W'_2, using its three-premise formulation.
    for x, y, r in product(range(len(s)), repeat=3):
        ddx = table[table[x][zero]][zero]
        ddy = table[table[y][zero]][zero]
        if (table[x][table[y][zero]] == one and
                table[table[ddx][x]][r] == one and
                table[table[ddy][y]][r] == one):
            assert r == one
    result, proper = direct_critical(table, one, zero)
    assert result and criterion(up, s)
    print(f'SOURCE_S n={n} points={len(up)} elements={len(s)} '
          f'cardinality_height_Wprime2=OK direct_critical={result} '
          f'proper_endomorphisms={len(proper)}', flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--max-points', type=int, default=6)
    parser.add_argument('--source-max-n', type=int, default=8)
    args = parser.parse_args()
    total = critical = normalizations = 0
    sizes = {}
    for n in range(args.max_points + 1):
        pc = fc = cc = nc = 0
        largest = 0
        for up in posets(n):
            pc += 1
            for s in families(up):
                table, one, zero = algebra(up, s)
                identities(up, s, table, one, zero)
                actual, proper = direct_critical(table, one, zero)
                predicted = criterion(up, s)
                assert actual == predicted, (up, s, actual, predicted)
                nc += normalization_checks(up, s, table, one, proper)
                cc += actual
                fc += 1
                largest = max(largest, len(s))
                sizes[len(s)] = sizes.get(len(s), 0) + 1
            if n == args.max_points and pc % 50 == 0:
                print(f'EXHAUSTIVE_PREFIX points={n} posets={pc} '
                      f'hereditary_families={fc} critical={cc}', flush=True)
        assert pc == (1, 1, 2, 5, 16, 63, 318)[n]
        total += fc
        critical += cc
        normalizations += nc
        print(f'EXHAUSTIVE points={n} unlabeled_posets={pc} '
              f'hereditary_families={fc} critical={cc} max_elements={largest} '
              f'proper_separating_d_normalizations={nc}', flush=True)
    for n in range(2, args.source_max_n + 1):
        source_checks(n)
    print(f'TOTAL representations={total} critical={critical} '
          f'fold_normalizations={normalizations}', flush=True)
    print('ALGEBRA_SIZE_HISTOGRAM ' + str(dict(sorted(sizes.items()))), flush=True)
    print('ALL_CHECKS_OK', flush=True)


if __name__ == '__main__':
    main()
