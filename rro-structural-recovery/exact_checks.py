"""Exact finite evidence for the literal RRO source/capacity/label evaluator.

Python >=3.12, standard library only. Run: python3 exact_checks.py --output FILE
Source definitions: UJB.7502 and UJB.7701; labels use MSB-first coordinates.
No full hidden injection is enumerated. Small-graph homomorphism inversion is
used as a second route for the <=4-sample evaluator, not as the general proof.
"""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import cache
from itertools import product, combinations, permutations
import json
from math import factorial
import platform
import sys
from pathlib import Path

assert sys.version_info >= (3, 12)
M = 1944
V = tuple(x for x in product(range(3), repeat=2) if any(x))
Z = tuple(product(range(3), repeat=2))
W = tuple(product(range(3), repeat=3))
LEFT = tuple(product(V, Z))
RIGHT = tuple(product(W, Z))


def bilinear(i, v, w):
    x, y = v
    a, b, c = w
    return ((x*a+y*b) % 3, (x*b+y*c) % 3) if i == 0 else ((x*a-y*b) % 3, (y*a+x*b) % 3)


def falling(n, k):
    return factorial(n)//factorial(n-k) if 0 <= k <= n else 0


def literal(i):
    return tuple(frozenset(j for j, (w, z) in enumerate(RIGHT) if bilinear(i, v, w) == s) for v, s in LEFT)


def twin_geometry(i, adj):
    lclasses = defaultdict(list)
    rclasses = defaultdict(list)
    for j, row in enumerate(adj):
        lclasses[row].append(j)
    cols = tuple(frozenset(a for a, row in enumerate(adj) if t in row) for t in range(243))
    for j, col in enumerate(cols):
        rclasses[col].append(j)
    lc, rc = tuple(lclasses.values()), tuple(rclasses.values())
    li = {x: j for j, C in enumerate(lc) for x in C}
    ri = {x: j for j, C in enumerate(rc) for x in C}
    edges = {(li[a], ri[t]) for a, row in enumerate(adj) for t in row}
    # Explicit geometric class keys, obtained independently from the neighborhoods.
    if i == 0:
        leftkeys = []
        for v, s in LEFT:
            pivot = next(x for x in v if x)
            inv = 1 if pivot == 1 else 2
            leftkeys.append(tuple(inv*x % 3 for x in v+s))
        rightkeys = [w for w, z in RIGHT]
        directions = set()
        for C in lc:
            points = {RIGHT[t][0] for t in adj[C[0]]}
            assert len(points) == 3
            a, b = tuple(points)[:2]
            delta = tuple((x-y) % 3 for x, y in zip(a, b))
            pivot = next(x for x in delta if x)
            directions.add(tuple(x*(1 if pivot == 1 else 2) % 3 for x in delta))
        assert len(directions) == 4
    else:
        leftkeys = []
        for (x, y), (a, b) in LEFT:
            # s/v = s*conjugate(v)/(x^2+y^2) in F_9.
            norm = (x*x+y*y) % 3
            inv = 1 if norm == 1 else 2
            leftkeys.append(((a*x+b*y)*inv % 3, (b*x-a*y)*inv % 3))
        rightkeys = [w[:2] for w, z in RIGHT]
        assert all(leftkeys[a] == rightkeys[t] for a, row in enumerate(adj) for t in row)
    for keys, classes in ((leftkeys, lc), (rightkeys, rc)):
        assert all(len({keys[x] for x in C}) == 1 for C in classes)
        assert len({keys[C[0]] for C in classes}) == len(classes)
    unseen = {('L', j) for j in range(len(lc))} | {('R', j) for j in range(len(rc))}
    components = []
    while unseen:
        todo, seen = [next(iter(unseen))], set()
        while todo:
            node = todo.pop()
            if node in seen:
                continue
            seen.add(node)
            side, j = node
            todo.extend(('R', b) for a, b in edges if side == 'L' and a == j)
            todo.extend(('L', a) for a, b in edges if side == 'R' and b == j)
        unseen -= seen
        components.append([sum(s == 'L' for s, _ in seen), sum(s == 'R' for s, _ in seen)])
    return lc, rc, edges, cols, {'left_classes': len(lc), 'left_sizes': dict(Counter(map(len, lc))), 'right_classes': len(rc), 'right_sizes': dict(Counter(map(len, rc))), 'quotient_edges': len(edges), 'components': sorted(components)}


def capacity_embeddings(r, s, edges, geometry):
    """Direct type assignment with remaining capacities; generic finite graph."""
    lc, rc, incidences, _, _ = geometry
    leftneighbors = [set() for _ in lc]
    rightneighbors = [set() for _ in rc]
    for a, b in incidences:
        leftneighbors[a].add(b)
        rightneighbors[b].add(a)
    graph = [set() for _ in range(r+s)]
    for a, b in edges:
        graph[a].add(r+b)
        graph[r+b].add(a)
    usedl, usedr = [0]*len(lc), [0]*len(rc)
    assigned = {}

    def domain(x):
        neighbors = rightneighbors if x < r else leftneighbors
        n = len(lc) if x < r else len(rc)
        result = set(range(n))
        for y in graph[x]:
            if y in assigned:
                result &= neighbors[assigned[y]]
        return result

    def go():
        if len(assigned) == r+s:
            return 1
        x = min((x for x in range(r+s) if x not in assigned), key=lambda x: (len(domain(x)), -len(graph[x])))
        used, classes = (usedl, lc) if x < r else (usedr, rc)
        total = 0
        for typ in domain(x):
            available = len(classes[typ])-used[typ]
            if available:
                assigned[x] = typ
                used[typ] += 1
                total += available*go()
                used[typ] -= 1
                del assigned[x]
        return total
    return go()


@cache
def partitions(n):
    if n == 0:
        return ((),)
    out = []
    for p in partitions(n-1):
        for a in range(max(p, default=-1)+2):
            out.append(p+(a,))
    return tuple(out)


def blocks(p):
    return tuple(tuple(j for j, x in enumerate(p) if x == a) for a in range(max(p, default=-1)+1))


def mobius(p):
    answer = 1
    for B in blocks(p):
        answer *= (-1)**(len(B)-1)*factorial(len(B)-1)
    return answer


@cache
def hom_at_most_four(i, r, s, edges):
    graph = [set() for _ in range(r+s)]
    for a, b in edges:
        graph[a].add(r+b)
        graph[r+b].add(a)
    unseen, total = set(range(r+s)), 1
    while unseen:
        todo, seen = [next(iter(unseen))], set()
        while todo:
            x = todo.pop()
            if x not in seen:
                seen.add(x)
                todo.extend(graph[x]-seen)
        unseen -= seen
        a, b = sum(x < r for x in seen), sum(x >= r for x in seen)
        e = sum(len(graph[x]) for x in seen)//2
        if e == 0:
            total *= 72 if a else 243
        elif e == len(seen)-1:
            total *= M*8**(a-1)*27**(b-1)
        else:
            assert (a, b, e) == (2, 2, 4)
            total *= (209952, 419904)[i]
    return total


@cache
def injection_by_inversion(i, r, s, edges):
    total = 0
    for p in partitions(r):
        for q in partitions(s):
            collapsed = tuple(sorted({(p[a], q[b]) for a, b in edges}))
            total += mobius(p)*mobius(q)*hom_at_most_four(i, max(p)+1, max(q)+1, collapsed)
    return total


def allowed(assignment, width):
    return tuple(v for v in range(1 << width) if all((v >> (width-1-h)) & 1 == bit for h, bit in assignment))


def merge_masks(B, cache, offset, width):
    result = {}
    for j in B:
        for h, bit in cache[j]:
            if offset <= h < offset+width:
                h -= offset
                if h in result and result[h] != bit:
                    return None
                result[h] = bit
    return tuple(sorted(result.items()))


@cache
def label_dp(masks, width):
    if any(m is None for m in masks):
        return 0
    r, full = len(masks), (1 << len(masks))-1
    sets = [set(allowed(m, width)) for m in masks]
    f = [0]*(full+1)
    f[0] = 1
    for label in range(1 << width):
        g = f.copy()
        for S in range(1, full+1):
            for b in range(r):
                if S & (1 << b) and label in sets[b]:
                    g[S] += f[S ^ (1 << b)]
        f = g
    return f[full]


@cache
def label_inclusion_exclusion(masks, width):
    """Independent count via collision partitions and intersections of masks."""
    if any(m is None for m in masks):
        return 0
    total = 0
    for p in partitions(len(masks)):
        term = mobius(p)
        for B in blocks(p):
            merged = {}
            conflict = False
            for b in B:
                for h, bit in masks[b]:
                    if h in merged and merged[h] != bit:
                        conflict = True
                    merged[h] = bit
            term *= 0 if conflict else 1 << (width-len(merged))
        total += term
    return total


def evidence(cache):
    """At most four samples; actual graphs, no full-injection enumeration."""
    n = len(cache)
    assert n <= 4
    values = [F(0), F(0)]
    for p in partitions(n):
        r = max(p)+1
        lm = tuple(merge_masks(B, cache, 0, 7) for B in blocks(p))
        l = label_dp(lm, 7)
        assert l == label_inclusion_exclusion(lm, 7)
        if not l:
            continue
        for q in partitions(n):
            s = max(q)+1
            rm = tuple(merge_masks(B, cache, 7, 8) for B in blocks(q))
            rr = label_dp(rm, 8)
            assert rr == label_inclusion_exclusion(rm, 8)
            if not rr:
                continue
            edges = tuple(sorted(set(zip(p, q))))
            for i in (0, 1):
                values[i] += F(injection_by_inversion(i, r, s, edges)*l*rr, falling(128, r)*falling(256, s)*M**n)
    return tuple(values)


def complete(a, t):
    return tuple((h, (a >> (6-h)) & 1) for h in range(7))+tuple((7+h, (t >> (7-h)) & 1) for h in range(8))


def conjoin(cache, j, h, bit):
    out = list(cache)
    old = dict(out[j])
    if h in old and old[h] != bit:
        return None
    old[h] = bit
    out[j] = tuple(sorted(old.items()))
    return tuple(out)


def mix(values, theta):
    return (1-theta)*values[0]+theta*values[1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    options = parser.parse_args()
    result = {'environment': {'python': sys.version, 'executable': sys.executable, 'platform': platform.platform()}, 'sources': []}
    geoms = []
    adjs = []
    for i in (0, 1):
        adj = literal(i)
        adjs.append(adj)
        geom = twin_geometry(i, adj)
        geoms.append(geom)
        hist = Counter(len(adj[a] & adj[b]) for a, b in combinations(range(72), 2))
        tri = Counter(len(adj[a] & adj[b] & adj[c]) for a, b, c in combinations(range(72), 3))
        rightpair = Counter(len(geom[3][u] & geom[3][v]) for u, v in combinations(range(243), 2))
        assert hist == ({0: 1872, 9: 648, 27: 36}, {0: 2304, 27: 252})[i]
        assert max(tri) == (9, 27)[i]
        direct_square = sum(2*n*falling(h, 2) for h, n in hist.items())
        right_square = sum(2*n*falling(h, 2) for h, n in rightpair.items())
        assert direct_square == right_square == (143856, 353808)[i]
        direct_k310 = sum(6*n*falling(h, 10) for h, n in tri.items())
        formula_k310 = 0 if i == 0 else 9*falling(8, 3)*falling(27, 10)
        assert direct_k310 == formula_k310 == (0, 92575499208192000)[i]
        for m in (2, 9, 10, 27, 28):
            direct = sum(6*n*falling(h, m) for h, n in tri.items())
            formula = (27*falling(8, 3)*falling(9, m) if m >= 2 else None) if i == 0 else 9*falling(8, 3)*falling(27, m)
            if i == 0 and m == 2:
                formula = 27*falling(8, 3)*falling(9, 2)
            assert direct == formula
        # Actual-source proof of forced right twins, checked over all distinct pairs.
        assert all(cols1 == cols2 for u, cols1 in enumerate(geom[3]) for cols2 in geom[3][u+1:] if len(cols1 & cols2) >= 3)
        path_direct = sum((len(adj[a])-1)*(len(geom[3][t])-1) for a, row in enumerate(adj) for t in row)
        assert path_direct == 353808
        square_capacity = capacity_embeddings(2, 2, ((0, 0), (0, 1), (1, 0), (1, 1)), geom)
        path_capacity = capacity_embeddings(2, 2, ((0, 0), (0, 1), (1, 0)), geom)
        assert square_capacity == direct_square and path_capacity == path_direct
        # Connected subgraphs exercise capacity recurrence against independent inversion.
        for r, s, edges in [(1, 1, ((0, 0),)), (1, 3, ((0, 0), (0, 1), (0, 2))), (3, 1, ((0, 0), (1, 0), (2, 0)))]:
            assert capacity_embeddings(r, s, edges, geom) == injection_by_inversion(i, r, s, edges)
        five = sum(2*n*falling(h, 2)*6 for h, n in rightpair.items())
        six = sum(2*n*falling(h, 3) for h, n in rightpair.items())
        assert (five, six) == ((863136, 653184), (2122848, 2122848))[i]
        result['sources'].append({'source': i, **geom[4], 'edges': sum(map(len, adj)), 'left_pair_histogram': dict(sorted(hist.items())), 'left_triple_histogram': dict(sorted(tri.items())), 'right_pair_histogram': dict(sorted(rightpair.items())), 'P4': path_direct, 'K22': direct_square, 'K310': direct_k310, 'five_edges': five, 'K32': six})
    # Positive forcing and sharp nonforcing witnesses in literal coordinates.
    for adj in adjs:
        aa = [LEFT.index((v, (0, 0))) for v in ((1, 0), (2, 0), (0, 1), (0, 2))]
        uu = RIGHT.index(((0, 0, 0), (0, 0)))
        vv = RIGHT.index(((0, 0, 0), (1, 0)))
        assert all(uu in adj[a] and vv in adj[a] for a in aa)
    a, b, d = [LEFT.index((v, (0, 0))) for v in ((1, 0), (2, 0), (0, 1))]
    u, v = [RIGHT.index((w, (0, 0))) for w in ((0, 0, 0), (0, 0, 1))]
    assert all(t in adjs[0][x] for x, t in ((a, u), (a, v), (b, u), (b, v), (d, u))) and v not in adjs[0][d]
    # Independent label counting on actual 128/256 alphabets and partial masks.
    mask_cases = 0
    for width in (7, 8):
        for r in (1, 2, 3, 4):
            for trial in range(20):
                masks = tuple(tuple((h, (trial+b+h) % 2) for h in range(width) if (trial+3*b+5*h) % 4 == 0) for b in range(r))
                assert label_dp(masks, width) == label_inclusion_exclusion(masks, width)
                mask_cases += 1
    empty = ((),)*4
    assert evidence(empty) == (1, 1)
    base = (complete(0, 0), complete(0, 1), complete(1, 0), ())
    zbase = evidence(base)
    assert zbase == (F(353808, M**3*falling(128, 2)*falling(256, 2)),)*2
    closed = base[:3]+(complete(1, 1),)
    zclosed = evidence(closed)
    predictions = {}
    for theta in (F(0), F(1, 3), F(1, 2), F(1)):
        pred = mix(zclosed, theta)/mix(zbase, theta)
        assert pred == F(37+54*theta, 91*M)
        predictions[str(theta)] = str(pred)
    # Every single coordinate splits a partial evidence probability exactly.
    prefix = empty
    transitions = 0
    target = complete(0, 0)
    for h, bit in target:
        zero, one = conjoin(prefix, 0, h, 0), conjoin(prefix, 0, h, 1)
        z, z0, z1 = evidence(prefix), evidence(zero), evidence(one)
        assert all(z0[i]+z1[i] == z[i] for i in (0, 1))
        prefix = (zero, one)[bit]
        assert conjoin(prefix, 0, h, bit) == prefix
        assert conjoin(prefix, 0, h, 1-bit) is None
        transitions += 1
    assert prefix[0] == complete(0, 0)
    assert evidence(prefix) == (F(1, 128*256),)*2
    # Regression 13: exact full-past versus coarsened-past first-bit laws.
    regressions = {}
    for bit, expected in ((0, F(575, 1143)), (1, F(568, 1143))):
        past = (complete(bit*64, 0), (), (), ())
        nxt = conjoin(past, 1, 0, 0)
        zp, zn = evidence(past), evidence(nxt)
        assert tuple(zn[i]/zp[i] for i in (0, 1)) == (expected,)*2
        regressions[str(bit)] = str(expected)
    # Actual old-slot permutations have identical likelihoods but different label factors.
    d0011 = tuple(tuple(sorted(((0, 0), (1, bit), (7, 0)))) for bit in (0, 0, 1, 1))
    d0101 = tuple(tuple(sorted(((0, 0), (1, bit), (7, 0)))) for bit in (0, 1, 0, 1))
    assert evidence(d0011) == evidence(d0101)
    # One fixed permitted full generator independently verifies the strict witness.
    old_slot = []
    for i, adj in enumerate(adjs):
        a0 = LEFT.index(((1, 0), (0, 0)))
        a1 = LEFT.index(((0, 1), (0, 0)))
        t0 = RIGHT.index(((0, 0, 0), (0, 0)))
        lm = {a0: 0, a1: 33}
        remaining = iter(x for x in range(128) if x not in (0, 32, 33))
        for a in range(72):
            if a not in lm:
                lm[a] = next(remaining)
        rm = {t0: 0}
        for t, lab in zip((t for t in range(243) if t != t0), range(1, 243)):
            rm[t] = lab
        counts, suffixes = [0, 0], [Counter(), Counter()]
        for a, row in enumerate(adj):
            for t in row:
                if lm[a] < 64 and rm[t] < 128:
                    bit = (lm[a] >> 5) & 1
                    counts[bit] += 1
                    suffixes[bit][(lm[a] % 32, rm[t])] += 1
        norm = sum((F(suffixes[0][x], counts[0])-F(suffixes[1][x], counts[1]))**2 for x in suffixes[0].keys() | suffixes[1].keys())
        assert counts == ([447, 449], [458, 438])[i]
        assert norm == (F(716, 200703), F(367, 100302))[i]
        old_slot.append({'source': i, 'mask_counts': counts, 'squared_suffix_difference': str(norm)})
    # Regression 11: actual uniform left-image bucket prior, independently summed.
    from math import comb
    bucket_z = bucket_pp = bucket_pr = F(0)
    bucket_vectors = 0
    total_images = comb(128, 72)
    for n00 in range(33):
        for n01 in range(33):
            for n10 in range(33):
                n11 = 72-n00-n01-n10
                if not 0 <= n11 <= 32:
                    continue
                bucket_vectors += 1
                weight = F(comb(32, n00)*comb(32, n01)*comb(32, n10)*comb(32, n11), total_images)
                pp = F(n00*n00*n11*n11, 72**4)
                rr = F(n01*n01*n10*n10, 72**4)
                bucket_z += weight*pp
                bucket_pp += weight*pp*pp
                bucket_pr += weight*pp*rr
    assert bucket_vectors == 22109
    assert bucket_z == F(1318541, 333298800)
    page_c = tuple(((0, bit), (1, bit)) for bit in (0, 0, 1, 1))
    page_cp = tuple(((0, bit), (1, 1-bit)) for bit in (0, 0, 1, 1))
    assert evidence(page_c) == evidence(page_cp) == (bucket_z, bucket_z)
    assert bucket_pp/bucket_z == F(59793742630939, 13730035262017536)
    assert bucket_pr/bucket_z == F(32023312367, 8938825040376)
    # Full stationary Read classification versus independent recursive term readout.
    alpha, beta = 'alpha', 'beta'
    C = (beta, alpha)
    gadgets = ((beta, C), (C, beta))
    marker = ((C, alpha), alpha)
    table = {(): ('b', 'b'), (0,): (beta, 'b'), (1,): ('b', beta),
             (0, 0): ('absent', beta), (0, 1): ('absent', alpha),
             (1, 0): (beta, 'absent'), (1, 1): (alpha, 'absent')}
    def tree_read(t, word):
        for direction in word:
            if not isinstance(t, tuple):
                return 'absent'
            t = t[direction]
        return 'b' if isinstance(t, tuple) else t
    def classified(bits, word):
        k = len(bits)
        for b in range(k):
            if b >= len(word):
                return 'b'
            if word[b] == 0:
                suffix = word[b+1:]
                return table[suffix][bits[b]] if suffix in table else 'absent'
        return tree_read(marker, word[k:])
    words = set()
    for n in range(11):
        words.update(product((0, 1), repeat=n))
    for b in range(60):
        for n in range(5):
            words.update((1,)*b+(0,)+u for u in product((0, 1), repeat=n))
    for n in range(8):
        words.update((1,)*60+u for u in product((0, 1), repeat=n))
    read_checks = 0
    for bits in ((0,)*60, (1,)*60, tuple(j % 2 for j in range(60)), tuple((j*j+j//3) % 2 for j in range(60))):
        tree = marker
        for bit in reversed(bits):
            tree = (gadgets[bit], tree)
        for word in words:
            assert tree_read(tree, word) == classified(bits, word)
            read_checks += 1
    from exact_evaluator import normalizer as direct_normalizer, partition_words
    for n in range(1, 7):
        assert tuple(partition_words(n, n)) == partitions(n)
    capacity_cases = []
    almost_closed = closed[:3]+(tuple(x for x in closed[3] if x[0] != 14),)
    for D in (closed, almost_closed):
        direct = direct_normalizer(D)
        independent = evidence(D)
        assert direct == independent
        capacity_cases.append({'mask_sizes': list(map(len, D)), 'capacity_Z': list(map(str, direct)), 'independent_Z': list(map(str, independent))})
    result['direct_capacity_evaluator'] = {'comparisons': capacity_cases, 'partition_generator_N_range': [1, 6]}
    result['full_native_menu'] = {'independent_term_vs_classifier_comparisons': read_checks, 'sample_bit_count': 60, 'covers': 'all words through length 10, every slot suffix through length 4, marker suffixes through length 7'}
    result['joint_page_regression'] = {'bucket_vectors': bucket_vectors, 'Z': str(bucket_z), 'future_given_c': str(bucket_pp/bucket_z), 'future_given_c_prime': str(bucket_pr/bucket_z)}
    result['evaluator'] = {'partition_pairs_at_N4': len(partitions(4))**2, 'independent_label_cases': mask_cases, 'partial_to_complete_updates': transitions, 'mixed_missing_pair': predictions, 'normalizer_P4': list(map(str, zbase)), 'normalizer_square': list(map(str, zclosed)), 'past_evidence_regression': regressions, 'old_slot_regression': old_slot, 'null_domain': 'Contradictory acquired bit conjunction is rejected; no conditional law is assigned.'}
    result['boundaries'] = ['Finite exact checks do not replace ordinary all-cache proofs.', 'No Lean compilation or full hidden injection enumeration.', 'No runtime-speed, polynomial-complexity, memory-optimality or world-novelty claim.']
    options.output.write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps({'exact_checks': 'passed', 'output': str(options.output), 'sources': 2, 'label_cases': mask_cases, 'completion_steps': transitions, 'mixed_prediction': predictions}, sort_keys=True))

if __name__ == '__main__':
    main()
