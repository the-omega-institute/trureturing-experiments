#!/usr/bin/env python3
"""Exact finite checks of Kostana Theorem 3.1 and the proposed pushout route.

Every group of order <=8 is included. Additional order-9 and order-16 cases
are listed below. No infinite-cardinal assertion follows from these checks.
"""
from functools import lru_cache
from itertools import product
from math import lcm


class Group:
    def __init__(self, moduli):
        self.moduli = tuple(moduli)
        self.elements = tuple(product(*(range(n) for n in moduli)))
        self.zero = (0,) * len(moduli)
        self.exp = lcm(*moduli) if moduli else 1

    def add(self, x, y):
        return tuple((a + b) % n for a, b, n in zip(x, y, self.moduli))

    def neg(self, x):
        return tuple(-a % n for a, n in zip(x, self.moduli))

    def mul(self, k, x):
        return tuple(k * a % n for a, n in zip(x, self.moduli))

    def generated(self, seeds):
        result = {self.zero}
        for a in seeds:
            result = {self.add(b, self.mul(k, a))
                      for b in result for k in range(self.exp)}
        return frozenset(result)

    @lru_cache(None)
    def subgroups(self):
        found = {frozenset({self.zero})}
        pending = list(found)
        while pending:
            subgroup = pending.pop()
            for a in self.elements:
                if a not in subgroup:
                    larger = self.generated(tuple(subgroup) + (a,))
                    if larger not in found:
                        found.add(larger)
                        pending.append(larger)
        return tuple(sorted(found, key=lambda h: (len(h), sorted(h))))

    @lru_cache(None)
    def multiples(self, k, subgroup=None):
        return frozenset(self.mul(k, a) for a in
                         (self.elements if subgroup is None else subgroup))

    @lru_cache(None)
    def pure(self, subgroup):
        # Multiplication maps repeat modulo the group exponent.
        return all(self.multiples(k, subgroup) ==
                   subgroup.intersection(self.multiples(k))
                   for k in range(1, self.exp + 1))


def homomorphisms(source, subgroup, target):
    """All maps, by generators and every equality h+k*a=h'+k'*a.

    No splitting hypothesis is used by the enumerator.
    """
    maps = [{source.zero: target.zero}]
    domain = frozenset({source.zero})
    for a in sorted(subgroup):
        if a in domain:
            continue
        enlarged = source.generated(tuple(domain) + (a,))
        next_maps = []
        for mapping in maps:
            for image in target.elements:
                extended = {}
                valid = True
                for h in domain:
                    for k in range(source.exp + 1):
                        x = source.add(h, source.mul(k, a))
                        y = target.add(mapping[h], target.mul(k, image))
                        if x in extended and extended[x] != y:
                            valid = False
                            break
                        extended[x] = y
                    if not valid:
                        break
                if valid:
                    assert set(extended) == set(enlarged)
                    next_maps.append(extended)
        maps, domain = next_maps, enlarged
    assert domain == subgroup
    assert all(mapping[source.add(x, y)] == target.add(mapping[x], mapping[y])
               for mapping in maps for x in subgroup for y in subgroup)
    return maps


def cosets(group, subgroup):
    return {x: frozenset(group.add(x, a) for a in subgroup)
            for x in group.elements}


def cocyclic_quotient(group, subgroup):
    size = len(group.elements) // len(subgroup)
    if size == 1:
        return True
    primes = [p for p in range(2, size + 1)
              if size % p == 0 and all(p % d for d in range(2, p))]
    q = cosets(group, subgroup)
    return len(primes) == 1 and any(
        len({q[group.mul(k, a)] for k in range(group.exp)}) == size
        for a in group.elements)


def fidelity(group):
    # Source Theorem 3.1: the maximal-quotient diagonal pure embedding.
    coordinates = []
    for a in group.elements:
        for k in range(group.exp + 1):
            candidates = [m for m in group.subgroups()
                          if a not in m and group.multiples(k) <= m]
            m = max(candidates, key=lambda h: (len(h), sorted(h))) if candidates \
                else frozenset(group.elements)
            assert cocyclic_quotient(group, m)
            coordinates.append(m)
    assert set.intersection(*(set(m) for m in coordinates)) == {group.zero}
    for k in range(1, group.exp + 1):
        divisible = [{group.add(b, h) for b in group.multiples(k) for h in m}
                     for m in coordinates]
        for a in group.elements:
            assert (a in group.multiples(k)) == all(a in q for q in divisible)
    # Source Lemma 1.2, for finite groups: pure subgroups are direct summands.
    for a in group.subgroups():
        if group.pure(a):
            assert any(len(a) * len(c) == len(group.elements) and
                       a.intersection(c) == {group.zero} for c in group.subgroups())
    return len(coordinates)


def check_diagram(b, a, d, f, all_h):
    ambient = Group(d.moduli + b.moduli)
    width = len(d.moduli)
    relation = frozenset(f[x] + b.neg(x) for x in a)
    assert ambient.generated(relation) == relation
    q = cosets(ambient, relation)
    i_d = {x: q[x + b.zero] for x in d.elements}
    i_b = {x: q[d.zero + x] for x in b.elements}
    assert len(set(i_d.values())) == len(d.elements)
    assert len(set(i_b.values())) == len(b.elements)
    assert all(i_d[f[x]] == i_b[x] for x in a)
    for k in range(1, ambient.exp + 1):
        kp = {q[ambient.mul(k, x)] for x in ambient.elements}
        for source, inclusion in [(d, i_d), (b, i_b)]:
            assert {inclusion[source.mul(k, x)] for x in source.elements} == \
                set(inclusion.values()).intersection(kp)
    points = set(q.values())
    assert len(points) * len(a) == len(d.elements) * len(b.elements)
    # A retraction P->D amounts to an extension h:B->D of f.
    h = next((h for h in all_h if all(h[x] == f[x] for x in a)), None)
    assert h is not None
    r = {}
    for x in ambient.elements:
        value = d.add(x[:width], h[x[width:]])
        assert q[x] not in r or r[q[x]] == value
        r[q[x]] = value
    assert all(r[i_d[x]] == x for x in d.elements)
    # Realize the fresh quotient K=B/A by cosets, independently of the pushout.
    quotient = cosets(b, a)
    extension = {x: (h[x], quotient[x]) for x in b.elements}
    target = set(product(d.elements, set(quotient.values())))
    assert len(set(extension.values())) == len(b.elements)
    assert all(extension[x] == (f[x], quotient[b.zero]) for x in a)
    for x in b.elements:
        for y in b.elements:
            z = b.add(x, y)
            sum_coset = quotient[b.add(next(iter(quotient[x])), next(iter(quotient[y])))]
            assert extension[z] == (d.add(h[x], h[y]), sum_coset)
    for k in range(1, ambient.exp + 1):
        k_target = {(d.mul(k, x), quotient[b.mul(k, next(iter(c)))]) for x, c in target}
        assert {extension[b.mul(k, x)] for x in b.elements} == \
            set(extension.values()).intersection(k_target)
    assert {(d.add(x, h[y]), quotient[y]) for x in d.elements for y in b.elements} == target
    return len(points)


def main():
    moduli = [(), (2,), (3,), (4,), (2, 2), (5,), (6,), (7,), (8,),
              (2, 4), (2, 2, 2), (9,), (3, 3), (16,), (2, 8), (4, 4),
              (2, 2, 4), (2, 2, 2, 2)]
    groups = [Group(m) for m in moduli]
    targets = [Group(m) for m in [(), (2,), (3,), (4,), (2, 2)]]
    coordinates = sum(fidelity(b) for b in groups)
    print('FIDELITY: Theorem 3.1 quotient diagonal and finite Lemma 1.2')
    print('B moduli:', moduli)
    print('groups:', len(groups), 'quotient coordinates:', coordinates)
    tested = points = 0
    for b in groups:
        for d in targets:
            all_h = homomorphisms(b, frozenset(b.elements), d)
            for a in b.subgroups():
                if len(a) > len(d.elements) or not b.pure(a):
                    continue
                for f in homomorphisms(b, a, d):
                    image = frozenset(f.values())
                    if len(image) == len(a) and d.pure(image):
                        points += check_diagram(b, a, d, f, all_h)
                        tested += 1
    c4, c2 = Group((4,)), Group((2,))
    assert not c4.pure(c4.generated([(2,)]))
    assert not any(h[(2,)] == (1,) for h in homomorphisms(c4, frozenset(c4.elements), c2))
    print('EXTENSION: every pure A<=B and pure embedding A->D in the stated ranges')
    print('D moduli:', [d.moduli for d in targets])
    print('pure diagrams:', tested, 'total quotient points:', points)
    print('NEGATIVE CONTROL: nonpure C2<=C4 identity into C2 has no homomorphic extension')
    print('LIMIT: finite checks do not prove regularity, infinite compactness, or back-and-forth')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
