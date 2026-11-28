"""Exact ML-A expressions, independent label recurrence and finite source checks.

These are arithmetic validation programs; their host execution is NOT charged as
ML-T/B execution. Cost counters here are not LAB-resource coordinates.
"""
from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from functools import cache
from itertools import combinations, product


def partitions(n):
    def go(prefix):
        if len(prefix) == n:
            yield prefix
        else:
            for b in range(max(prefix, default=-1)+2):
                yield from go(prefix+(b,))
    yield from go(())


@cache
def labels(masks, width):
    if any(mask is None for mask in masks):
        return 0
    size = 1 << len(masks)
    old = [1]+[0]*(size-1)
    for label in range(1 << width):
        new = old.copy()
        for subset in range(1, size):
            for b, restrictions in enumerate(masks):
                if subset & (1 << b) and all((label >> (width-1-h)) & 1 == value for h, value in restrictions):
                    new[subset] += old[subset ^ (1 << b)]
        old = new
    return old[-1]


def conjoin(a, b):
    if a is None or b is None:
        return None
    assert len(a) == len(b)
    out = []
    for x, y in zip(a, b):
        d = dict(x)
        for h, value in y:
            if h in d and d[h] != value:
                return None
            d[h] = value
        out.append(tuple(sorted(d.items())))
    return tuple(out)


class Evaluator:
    def __init__(self, supplier, general):
        self.old, self.general = supplier, general
        self.normalizers = {}
        self.counts = Counter()

    def normalizer(self, data):
        if data is None:
            return F(0), F(0)
        if data in self.normalizers:
            return self.normalizers[data]
        n, old = len(data), self.old
        assert n >= 4
        # Unread original indices integrate to one under the iid edge draws
        # conditional on the ONE shared injection; they are not new injections.
        active = tuple(row for row in data if row)
        if len(active) <= 4:
            result = self.small_formula(data)
            if active:
                assert result == old.evidence(active)
            else:
                assert result == (F(1), F(1))
        elif all(row == data[0] and len(row) == 15 for row in data):
            result = (F(1, 1944**(n-1)*128*256),)*2
        elif len({tuple((h, v) for h, v in row if h < 7) for row in data
                  if sum(h < 7 for h, _ in row) == 7}) > 72:
            result = F(0), F(0)
        else:
            # Retain unrestricted mathematical domain; never silently truncate N.
            result = self.general.normalizer(data)
        self.normalizers[data] = result
        self.counts['normalizers'] += 1
        return result

    def small_formula(self, original):
        old, total = self.old, [F(0), F(0)]
        # This validation route is independent of ML-A's enumeration cost. At
        # <=4 informative indices use the mature <=4-edge inversion supplier.
        data = tuple(row for row in original if row)
        if not data:
            return F(1), F(1)
        n = len(data)
        for p in partitions(n):
            r = max(p)+1
            left = labels(tuple(old.merge_masks(b, data, 0, 7) for b in old.blocks(p)), 7)
            for q in partitions(n):
                self.counts['partition_pairs'] += 1
                s = max(q)+1
                right = labels(tuple(old.merge_masks(b, data, 7, 8) for b in old.blocks(q)), 8)
                if not left or not right:
                    self.counts['zero_pruned_pairs'] += 1
                    continue
                edges = tuple(sorted(set(zip(p, q))))
                for i in (0, 1):
                    j = old.injection_by_inversion(i, r, s, edges)
                    total[i] += F(j*left*right, 1944**n*old.falling(128, r)*old.falling(256, s))
        return tuple(total)

    def evaluate(self, n, theta, data, cylinders):
        assert n == len(data) and n >= 4 and isinstance(theta, F) and 0 <= theta <= 1
        z0, z1 = self.normalizer(data)
        z = (1-theta)*z0+theta*z1
        values, joint = [z0, z1, z], []
        for c in cylinders:
            a, b = self.normalizer(conjoin(data, c))
            m = (1-theta)*a+theta*b
            values.extend((a, b, m))
            joint.append(m)
        return values, None if z == 0 else [m/z for m in joint]


def repeated_full_direct(old, data):
    """Actual ML-A partition and zero-pruning schedule at N=4.

    The surviving pair has r=s=1; enumerate *all* type maps, not a
    precomputed answer. No compiled tape execution is asserted here.
    """
    assert len(data) == 4 and all(row == data[0] and len(row) == 15 for row in data)
    total, counts = [F(0), F(0)], Counter()
    geometry = [old.twin_geometry(i, old.literal(i)) for i in (0, 1)]
    for p in partitions(4):
        r = max(p)+1
        l = labels(tuple(old.merge_masks(b, data, 0, 7) for b in old.blocks(p)), 7)
        for q in partitions(4):
            counts['partition_pairs'] += 1
            s = max(q)+1
            right = labels(tuple(old.merge_masks(b, data, 7, 8) for b in old.blocks(q)), 8)
            if not l or not right:
                counts['zero_pruned_pairs'] += 1
                continue
            assert r == s == 1
            edges = tuple(zip(p, q))  # Keep every original occurrence.
            for i, (lc, rc, incidences, _, _) in enumerate(geometry):
                j = 0
                for lt in product(range(len(lc)), repeat=r):
                    for rt in product(range(len(rc)), repeat=s):
                        counts[f'G{i}_type_maps'] += 1
                        weight = 1
                        usedl, usedr = Counter(), Counter()
                        for t in lt:
                            weight *= max(0, len(lc[t])-usedl[t])
                            usedl[t] += 1
                        for t in rt:
                            weight *= max(0, len(rc[t])-usedr[t])
                            usedr[t] += 1
                        for a, b in edges:
                            counts[f'G{i}_incidence_tests'] += 1
                            if (lt[a], rt[b]) not in incidences:
                                weight = 0
                        j += weight
                assert j == 1944
                total[i] += F(j*l*right, 1944**4*old.falling(128, r)*old.falling(256, s))
    return tuple(total), dict(counts)


def capacity_witness(old):
    """Literal common-neighbor counting, independent of type-map formula."""
    injections, triples = [], []
    for i in (0, 1):
        rows = old.literal(i)
        j, positive, histogram = 0, 0, Counter()
        for a, b, c in combinations(range(72), 3):
            size = len(rows[a] & rows[b] & rows[c])
            histogram[size] += 1
            j += 6*old.falling(size, 10)
            positive += int(size >= 10)
        injections.append(j)
        triples.append({'triples': sum(histogram.values()), 'positive_triples': positive,
                        'common_neighbors': dict(sorted(histogram.items()))})
    assert injections[0] == 0
    assert injections[1] == 9*old.falling(8, 3)*old.falling(27, 10)
    z1 = F(injections[1], 1944**30*old.falling(128, 3)*old.falling(256, 10))
    return (F(0), z1), {'literal_lift_counts': injections, 'triple_census': triples,
                        'Z0': '0', 'Z1': str(z1), 'N': 30}
