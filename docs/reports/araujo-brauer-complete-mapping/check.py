#!/usr/bin/env python3
"""Exact finite checks for the existence clause of arXiv:2608.25092, Problem 15.11.

This checks the proposed general argument's certificates, not its Lean proof.
No numerical approximation, SAT oracle, or third-party Python dependency is used.
"""
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations
from math import comb, factorial


class UnionFind:
    def __init__(self, n):
        self.p = list(range(n))

    def root(self, x):
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def join(self, x, y):
        self.p[self.root(x)] = self.root(y)


@lru_cache(None)
def matchings(vertices, partial=False):
    """A matching is a tuple of ordered pairs; unmatched vertices are singletons."""
    if not vertices:
        return ((),)
    x, *tail = vertices
    answer = []
    if partial:
        answer.extend(matchings(tuple(tail), True))
    for i, y in enumerate(tail):
        for rest in matchings(tuple(tail[:i] + tail[i + 1:]), partial):
            answer.append(((x, y),) + rest)
    return tuple(answer)


def telephone(k):
    a, b = 1, 1
    for j in range(2, k + 1):
        a, b = b, b + (j - 1) * a
    return 1 if k == 0 else b


def pairing_count(k):
    assert k >= 0 and k % 2 == 0
    ans = 1
    for j in range(1, k, 2):
        ans *= j
    return ans


def half_diagrams(n, r, partial=False):
    out = []
    for through in combinations(range(1, n + 1), r):
        complement = tuple(x for x in range(1, n + 1) if x not in through)
        out.extend((through, p) for p in matchings(complement, partial))
    return out


def transport(n, a, b):
    """Sandwich entry: ordered through-vertex transport; None means lost rank."""
    uf = UnionFind(n + 1)
    for x, y in a[1] + b[1]:
        uf.join(x, y)
    targets = {}
    for j, y in enumerate(b[0]):
        targets.setdefault(uf.root(y), []).append(j)
    p = []
    for x in a[0]:
        js = targets.get(uf.root(x), [])
        if len(js) != 1:
            return None
        p.append(js[0])
    return tuple(p) if len(set(p)) == len(p) else None


def sign(p):
    assert p is not None
    return (-1) ** sum(p[i] > p[j] for i in range(len(p)) for j in range(i + 1, len(p)))


def complete_mappings(table, count_all=False):
    """Exact cover of row, column and product constraints; returns witnesses."""
    n = len(table)
    answers = []
    alpha = [-1] * n

    def visit(rows, used_columns, used_products):
        if not rows:
            witness = tuple(alpha)
            assert sorted(witness) == list(range(n))
            assert sorted(table[x][witness[x]] for x in range(n)) == list(range(n))
            answers.append(witness)
            return not count_all
        best, options = None, None
        for x in sorted(rows):
            opts = [(y, table[x][y]) for y in range(n)
                    if y not in used_columns and table[x][y] not in used_products]
            if options is None or len(opts) < len(options):
                best, options = x, opts
            if not opts:
                return False
        for y, z in options:
            alpha[best] = y
            if visit(rows - {best}, used_columns | {y}, used_products | {z}):
                return True
        alpha[best] = -1
        return False

    visit(set(range(n)), set(), set())
    return answers


def check_source_fidelity():
    # Source Section 6, proposition labelled prop:brandt-latin:
    # all complete mappings of the matrix-unit semigroup correspond to Latin squares.
    # (The exact label/number is documented in README.)
    for k, expected in [(2, 2), (3, 12)]:
        elems = [None] + list((i, j) for i in range(k) for j in range(k))
        index = {x: i for i, x in enumerate(elems)}
        table = [[index[(x[0], y[1])] if x is not None and y is not None
                  and x[1] == y[0] else 0 for y in elems] for x in elems]
        found = complete_mappings(table, True)
        assert len(found) == expected
        for a in found:
            assert a[0] == 0
            square = [[elems[a[index[(i, j)]]][1] for j in range(k)] for i in range(k)]
            assert all(elems[a[index[(i, j)]]][0] == j for i in range(k) for j in range(k))
            assert all(sorted(row) == list(range(k)) for row in square)
            assert all(sorted(square[i][j] for i in range(k)) == list(range(k)) for j in range(k))
        print(f'FIDELITY matrix-unit order={k}: complete mappings={len(found)}=Latin squares')


BASE2 = [((2,4),(1,4)), ((3,4),(3,5)), ((1,5),(1,6)),
         ((4,5),(2,5)), ((4,6),(2,6)), ((5,6),(3,6))]
BASE3 = [((1,2,3),(1,2,5)), ((1,4,5),(1,3,5)), ((2,3,5),(2,3,6)),
         ((2,4,5),(2,4,6)), ((3,4,5),(3,4,6)), ((3,5,6),(1,3,6)),
         ((4,5,6),(1,4,6)), ((1,2,6),(1,2,7)), ((1,5,6),(1,5,7)),
         ((2,5,6),(2,5,7)), ((2,6,7),(1,6,7)), ((5,6,7),(3,5,7)),
         ((3,6,7),(1,3,7)), ((2,3,7),(2,4,7)), ((4,5,7),(1,4,7)),
         ((4,6,7),(3,4,7))]
# Each column is a partition of the j-subsets of four new vertices into adjacent pairs.
NEW_PAIRS = {
    1: [((0,), (1,)), ((2,), (3,))],
    2: [((0,1),(0,2)), ((0,3),(1,3)), ((1,2),(2,3))],
    3: [((0,1,2),(0,1,3)), ((0,2,3),(1,2,3))],
}


def triangle(r):
    common = () if r == 2 else (4,)
    return [tuple(sorted(t + common)) for t in [(1,2),(1,3),(2,3)]]


def through_pairs(n, r):
    base = 6 if r == 2 else 7
    assert n >= base and (n - base) % 4 == 0
    pairs = list(BASE2 if r == 2 else BASE3)
    for old in range(base, n, 4):
        for j in range(1, r + 1):
            for a in combinations(range(1, old + 1), r - j):
                for b, c in NEW_PAIRS[j]:
                    pairs.append((a + tuple(old + 1 + x for x in b),
                                  a + tuple(old + 1 + x for x in c)))
    return pairs


def check_subset_pairings():
    for j, pairs in NEW_PAIRS.items():
        flat = [a for pair in pairs for a in pair]
        assert sorted(flat) == list(combinations(range(4), j))
        assert all(len(set(a) & set(b)) == j - 1 for a, b in pairs)
    for r, base in [(2, 6), (3, 7)]:
        for n in range(base, 40, 4):
            pairs = through_pairs(n, r)
            flat = [a for pair in pairs for a in pair]
            assert len(set(flat)) == len(flat)
            assert set(flat) | set(triangle(r)) == set(combinations(range(1, n+1), r))
            assert not set(flat) & set(triangle(r))
            assert all(len(set(a) & set(b)) == r - 1 for a, b in pairs)
        print(f'SUBSETS rank={r}: base={base}, induction checked through n={n}')


def check_counts_and_obstructions():
    for k in range(11):
        assert len(matchings(tuple(range(k)), True)) == telephone(k)
        if k >= 2:
            assert telephone(k) % 2 == 0
        if k % 2 == 0:
            assert len(matchings(tuple(range(k)))) == pairing_count(k)
    for n in range(1, 41):
        for r in [2, 3]:
            if r < n:
                m = comb(n, r) * telephone(n-r)
                assert n < 4 or m % 2 == 0
            if r < n and (n-r) % 2 == 0:
                m = comb(n, r) * pairing_count(n-r)
                expected_odd = (r == 2 and n % 4 == 2) or (r == 3 and n % 4 == 3)
                assert bool(m % 2) == expected_odd
    for n in [2, 3]:
        group = list(permutations(range(n)))
        index = {p:i for i,p in enumerate(group)}
        table = [[index[tuple(q[p[i]] for i in range(n))] for q in group] for p in group]
        assert not complete_mappings(table)
    print('COUNTS telephone k=0..10; index parity n=1..40; S2/S3 obstructions: OK')


def check_small_supports():
    for partial in [False, True]:
        for n in range(1, 7):
            ranks = range(n + 1) if partial else range(n % 2, n+1, 2)
            for r in ranks:
                hs = half_diagrams(n, r, partial)
                expected = comb(n, r) * (telephone(n-r) if partial else pairing_count(n-r))
                assert len(hs) == expected
                # Full supports at n<=5; at n=6 sample actual diagrams for local compatibility.
                if n > 5:
                    sample = hs[:80]
                    for a in sample:
                        for b in sample:
                            if a[0] == b[0] or (not partial and len(set(a[0]) & set(b[0])) >= r-1):
                                assert transport(n, a, b) is not None
                    continue
                m = len(hs)
                support = [[transport(n, a, b) for b in hs] for a in hs]
                degrees = [sum(p is not None for p in row) for row in support]
                for i, a in enumerate(hs):
                    assert support[i][i] == tuple(range(r))
                    for j, b in enumerate(hs):
                        p = support[i][j]
                        q = support[j][i]
                        assert (p is None) == (q is None)
                        if p is not None:
                            assert all(q[p[k]] == k for k in range(r))
                        if a[0] == b[0] or (not partial and len(set(a[0]) & set(b[0])) >= r-1):
                            assert p is not None
                # The source's diagonal-slack construction needs only symmetry and loops.
                weights = [[Fraction(m-degrees[i]+1 if i == j else 1)
                            if support[i][j] is not None else Fraction(0)
                            for j in range(m)] for i in range(m)]
                assert all(sum(row) == m for row in weights)
                assert all(sum(weights[i][j] for i in range(m)) == m for j in range(m))
                assert all(weights[i][j] >= 1 for i in range(m) for j in range(m)
                           if support[i][j] is not None)
        print(f'SUPPORT {"PB" if partial else "B"}: full n=1..5, sampled n=6: OK')


def check_anchored_cycles():
    # Adjacent-through-set compatibility is specific to B, not PB.
    assert transport(2, ((1,), ()), ((2,), ())) is None
    print("BOUNDARY PB2 rank=1: adjacent through sets with only singleton complements are incompatible")
    for n, r in [(6,2), (7,3)]:
        hs = half_diagrams(n, r)
        fibres = {}
        for h in hs:
            fibres.setdefault(h[0], []).append(h)
        rho = {}
        def pair(a,b):
            assert a not in rho and b not in rho and a != b
            assert transport(n,a,b) is not None
            rho[a],rho[b] = b,a
        for a,b in through_pairs(n,r):
            for x,y in zip(fibres[a],fibres[b],strict=True):
                pair(x,y)
        ts = triangle(r)
        triple = [fibres[t][0] for t in ts]
        for i in range(3):
            rho[triple[i]] = triple[(i+1)%3]
        for t in ts:
            rest = fibres[t][1:]
            assert len(rest)%2 == 0
            for i in range(0,len(rest),2):
                pair(rest[i],rest[i+1])
        assert set(rho) == set(hs) == set(rho.values())
        triangle_transports = [transport(n,triple[i],triple[(i+1)%3]) for i in range(3)]
        assert triangle_transports == [tuple(range(r)),tuple(range(r)),(1,0)+tuple(range(2,r))]
        seen = set(); cycle_counts = Counter()
        for a in hs:
            if a in seen:
                continue
            x, length, product_sign = a,0,1
            while x not in seen:
                seen.add(x); length += 1
                product_sign *= sign(transport(n,x,rho[x]))
                x = rho[x]
            assert x == a
            assert product_sign == (-1)**length  # anchored equation, sigma=id, diagonal=id
            cycle_counts[length] += 1
        assert cycle_counts[3] == 1 and set(cycle_counts) <= {2,3}
        print(f'ANCHORED B{n} rank={r}: indices={len(hs)}, cycles={dict(sorted(cycle_counts.items()))}, triangle={triangle_transports}')


def check_n1_diagrams():
    # B1 is a singleton; PB1 consists of the identity and the all-singleton zero.
    for partial in [False,True]:
        ds = matchings((0,1),partial)
        def multiply(a,b):
            uf=UnionFind(3)
            for x,y in a:
                uf.join(x,y)
            for x,y in b:
                uf.join(x+1,y+1)
            return ((0,1),) if uf.root(0)==uf.root(2) else ()
        table = [[ds.index(multiply(a,b)) for b in ds] for a in ds]
        assert all(table[i][i] == i for i in range(len(ds)))
        assert complete_mappings(table)
    print('DEGREE1 B1/PB1: literal diagram products admit identity complete mappings')



def diagram_product(n, a, b):
    uf = UnionFind(3*n)
    for x,y in a:
        uf.join(x,y)
    for x,y in b:
        uf.join(x+n,y+n)
    exterior = list(range(n)) + list(range(2*n,3*n))
    groups = {}
    for i,x in enumerate(exterior):
        groups.setdefault(uf.root(x), []).append(i)
    assert all(len(xs) <= 2 for xs in groups.values())
    return tuple(sorted(tuple(xs) for xs in groups.values() if len(xs)==2))


def group_product(p,q):
    return tuple(q[p[i]] for i in range(len(p)))


def diagram_from_coordinates(n, a, g, b):
    pairs = [(x-1,y-1) for x,y in a[1]]
    pairs += [(x+n-1,y+n-1) for x,y in b[1]]
    pairs += [(x-1,b[0][g[k]]+n-1) for k,x in enumerate(a[0])]
    return tuple(sorted(pairs))


def check_full_small_monoids():
    group_maps = {}
    for r in [0,1,4]:
        group = list(permutations(range(r)))
        idx = {p:i for i,p in enumerate(group)}
        table = [[idx[group_product(p,q)] for q in group] for p in group]
        witnesses = complete_mappings(table)
        assert witnesses
        group_maps[r] = {p:group[witnesses[0][i]] for i,p in enumerate(group)}
    for partial in [False,True]:
        for n in range(1,5):
            diagrams = matchings(tuple(range(2*n)),partial)
            if n in [2,3]:
                # The unit ideal test, plus the exact S2/S3 search above, obstructs a mapping.
                def is_unit(d):
                    return len(d)==n and all(x<n<=y for x,y in d)
                units = {d for d in diagrams if is_unit(d)}
                assert len(units)==factorial(n)
                for a in diagrams:
                    for b in diagrams:
                        assert is_unit(diagram_product(n,a,b)) == (a in units and b in units)
                continue
            alpha = {}
            ranks = range(n+1) if partial else range(n%2,n+1,2)
            for r in ranks:
                hs = half_diagrams(n,r,partial)
                m = len(hs)
                group = list(permutations(range(r)))
                local = {}
                if r in [2,3]:
                    assert m%2==0
                    # Odd normal subgroup A_r: local complete mappings between its cosets.
                    for parity_in in [0,1]:
                        domain = [g for g in group if sign(g)==(-1)**parity_in]
                        for parity_out in [0,1]:
                            target = [g for g in group if sign(g)==(-1)**parity_out]
                            for values in permutations(target):
                                if len({group_product(g,f) for g,f in zip(domain,values)})==len(domain):
                                    local[parity_in,parity_out]=dict(zip(domain,values))
                                    break
                            assert (parity_in,parity_out) in local
                for i,a in enumerate(hs):
                    for lam,b in enumerate(hs):
                        for g in group:
                            source = diagram_from_coordinates(n,a,g,b)
                            if r in [2,3]:
                                eps = int(sign(g)==-1)
                                h = i//(m//2)
                                f = local[eps,h][g]
                                mu = (lam+2*i+eps)%m
                            else:
                                f = group_maps[r][g]
                                mu = (i+lam)%m
                            target = diagram_from_coordinates(n,b,f,hs[mu])
                            assert source not in alpha
                            alpha[source]=target
            assert set(alpha)==set(diagrams)==set(alpha.values())
            products = [diagram_product(n,x,alpha[x]) for x in diagrams]
            assert set(products)==set(diagrams) and len(set(products))==len(diagrams)
            print(f'FULL {"PB" if partial else "B"}{n}: explicit complete mapping on {len(diagrams)} diagrams')
        print(f'FULL {"PB" if partial else "B"}2/{"PB" if partial else "B"}3: unit-ideal obstruction verified on every diagram product')

def main():
    check_source_fidelity()
    check_n1_diagrams()
    check_counts_and_obstructions()
    check_subset_pairings()
    check_small_supports()
    check_anchored_cycles()
    check_full_small_monoids()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
