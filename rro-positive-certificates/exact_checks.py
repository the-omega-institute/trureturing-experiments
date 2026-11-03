"""Exact checks of the literal positive certificates; Python >= 3.12, stdlib.

Run from the experiment repository: python3 rro-positive-certificates/exact_checks.py
--output rro-positive-certificates/exact-results.json. Existing source helpers are
loaded from pinned Git objects in memory, so sparse rules need not be altered.
Uniform soundness is the ordinary proof in the accompanying public theory;
these finite checks do not enumerate all full injections or certify that proof.
"""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
from fractions import Fraction
import hashlib
from itertools import combinations, product
import json
from math import factorial
from pathlib import Path
import subprocess
import sys
import types

assert sys.version_info >= (3, 12)
BASE = '865f38257ead1b4e857ff915c23b0b04276bd7f3'
ROOT = Path(__file__).resolve().parent.parent
blob = subprocess.check_output(['git', '-C', str(ROOT), 'show',
                               BASE+':rro-structural-recovery/exact_checks.py'])
source = types.ModuleType('rro_source_checks')
exec(compile(blob, 'pinned:rro-structural-recovery/exact_checks.py', 'exec'), source.__dict__)
LEFT, RIGHT, M = source.LEFT, source.RIGHT, source.M


def literal_edges(i):
    return {(a, t) for a, row in enumerate(source.literal(i)) for t in row}


def index(D):
    rows, cols = defaultdict(set), defaultdict(set)
    for a, t in D:
        rows[a].add(t)
        cols[t].add(a)
    return rows, cols


def transport(E, lm, rm):
    return {(lm[a], rm[t]) for a, t in E}


def g0_witness(E):
    D, blocks = set(), [set(), set(), set()]
    for t, (w, z) in enumerate(RIGHT):
        if z == (0, 0):
            blocks[0].update((a, t) for a in range(72) if (a, t) in E)
        else:
            # v=(0,1) has kernel F(1,0,0); v=(1,0) has kernel F(0,0,1).
            vs = ((0, 1), (1, 0)) if w[0] in (0, 1) else ((0, 1),)
            for v in vs:
                s = source.bilinear(0, v, w)
                a = LEFT.index((v, s))
                blocks[1 if w[0] in (0, 1) else 2].add((a, t))
    D = set.union(*blocks)
    assert list(map(len, blocks)) == [216, 288, 72]
    assert len(D) == 576 and D <= E
    return D


def g0_complete(D):
    """Sufficient label-only completion on the promised positive domain."""
    rows, cols = index(D)
    if len(rows) != 72 or len(cols) != 243:
        return None
    anchors = sorted(t for t, ns in cols.items() if len(ns) == 8)
    if len(anchors) != 27 or len({frozenset(cols[t]) for t in anchors}) != 27:
        return None
    assigned = {t: t for t in anchors}
    counts = Counter(assigned.values())
    domains = {t: {u for u in anchors if ns <= cols[u]}
               for t, ns in cols.items() if t not in assigned}
    rounds = []
    while domains:
        singles = {t: next(iter(ds)) for t, ds in domains.items() if len(ds) == 1}
        if not singles:
            return None
        rounds.append(len(singles))
        for t, u in singles.items():
            assigned[t] = u
            counts[u] += 1
            if counts[u] > 9:
                return None
            del domains[t]
        full = {u for u in anchors if counts[u] == 9}
        domains = {t: ds-full for t, ds in domains.items()}
        if any(not ds for ds in domains.values()):
            return None
    if set(counts.values()) != {9}:
        return None
    return {(a, t) for t, u in assigned.items() for a in cols[u]}, rounds


def components(D):
    rows, cols = index(D)
    unseen = {(0, a) for a in rows} | {(1, t) for t in cols}
    out = []
    while unseen:
        todo, seen = [min(unseen)], set()
        while todo:
            side, x = todo.pop()
            if (side, x) in seen:
                continue
            seen.add((side, x))
            todo.extend((1-side, y) for y in (rows[x] if side == 0 else cols[x]))
        unseen -= seen
        out.append((sorted(x for s, x in seen if s == 0),
                    sorted(x for s, x in seen if s == 1)))
    return out


def g1_witness(E):
    comps = components(E)
    assert len(comps) == 9
    D = set()
    for j, (aa, tt) in enumerate(comps):
        if j == 0:
            D.update(product(aa[:3], tt[:10]))
            D.update((a, tt[0]) for a in aa[3:])
            D.update((aa[0], t) for t in tt[10:])
        else:
            D.update((aa[0], t) for t in tt)
            D.update((a, tt[0]) for a in aa[1:])
    assert len(D) == 324 and D <= E
    return D


def g1_complete(D, both=True):
    cc = components(D)
    if len(cc) != 9 or any((len(a), len(t)) != (8, 27) for a, t in cc):
        return None
    if both:
        rows, _ = index(D)
        if not any(len(rows[a] & rows[b] & rows[c]) >= 10
                   for aa, _ in cc for a, b, c in combinations(aa, 3)):
            return None
    return {(a, t) for aa, tt in cc for a in aa for t in tt}


def swap(n, a, b):
    p = list(range(n))
    p[a], p[b] = p[b], p[a]
    return p


def occupancy(q, n):
    counts = [1]+[0]*q
    for m in range(n):
        counts = [(M-q+k)*counts[k]+(q-k+1)*(counts[k-1] if k else 0)
                  for k in range(q+1)]
    return Fraction(counts[q], M**n)


ALPHA, BETA, BRANCH, ABSENT = 'alpha', 'beta', 'branch', 'absent'
C = (BETA, ALPHA)
GADGETS = ((BETA, C), (C, BETA))
MARKER = ((C, ALPHA), ALPHA)


def tree(bits):
    t = MARKER
    for b in reversed(bits):
        t = (GADGETS[b], t)
    return t


def read(t, word):
    for d in word:
        if not isinstance(t, tuple):
            return ABSENT
        t = t[d]
    return BRANCH if isinstance(t, tuple) else t


def bits_of(pair):
    a, t = pair
    return tuple((a >> h) & 1 for h in range(6, -1, -1))+tuple((t >> h) & 1 for h in range(7, -1, -1))


def conditioned(E, mask):
    return {y for y in E if all(bits_of(y)[h] == b for h, b in mask.items())}


def native_checks(E, witness):
    prepared = sorted(witness)+[min(E), max(E)]
    q, N = len(witness), len(prepared)
    bits = tuple(b for y in prepared for b in bits_of(y))
    term = tree(bits)
    acquired, total_length, maximum = [], 0, 0
    for b in range(15*q):
        word = (1,)*b+(0, 0)
        reply = read(term, word)
        assert reply in (BETA, BRANCH)
        acquired.append(int(reply == BRANCH))
        total_length += len(word)
        maximum = max(maximum, len(word))
    decoded = {(sum(acquired[15*j+h] << (6-h) for h in range(7)),
                sum(acquired[15*j+7+h] << (7-h) for h in range(8))) for j in range(q)}
    assert decoded == witness
    assert total_length == (15*q)*(15*q+3)//2 and maximum == 15*q+1
    # Original q-th slot is partial; (q+1)-th slot stays wholly unread.
    partial, unread = q, q+1
    actual = prepared[partial]
    masks = [dict(enumerate(bits_of(y))) for y in prepared[:q]]+[{}, {}]
    for h in range(7):
        b = 15*partial+h
        masks[partial][h] = int(read(term, (1,)*b+(0, 0)) == BRANCH)
    support_partial = conditioned(E, masks[partial])
    assert len(support_partial) == 27 and actual in support_partial
    assert len(conditioned(E, masks[unread])) == 1944
    assert conditioned(E, masks[0]) == {prepared[0]}
    # Independent local gadget lookup versus full original-root tree on old,
    # partial and unread indices; constants, informative absences and deep absences.
    probes = 0
    for j in (0, partial, unread):
        for h in range(15):
            bit = bits[15*j+h]
            for length in range(4):
                for suffix in product((0, 1), repeat=length):
                    word = (1,)*(15*j+h)+(0,)+suffix
                    assert read(term, word) == read(GADGETS[bit], suffix)
                    probes += 1
    # Conditional masses over every actually possible partial-slot edge.
    updates = []
    for h in range(7, 15):
        before = conditioned(E, masks[partial])
        zero = conditioned(E, masks[partial] | {h: 0})
        one = conditioned(E, masks[partial] | {h: 1})
        assert zero.isdisjoint(one) and zero | one == before
        zprob = Fraction(len(zero), len(before))
        assert zprob+Fraction(len(one), len(before)) == 1
        updates.append(str(zprob))
        masks[partial][h] = bits_of(actual)[h]
    assert conditioned(E, masks[partial]) == {actual}
    # All masks remain at original indices, including unread-slot replay domains.
    assert masks[unread] == {} and len(masks[partial]) == 15
    assert source.conjoin(tuple(tuple(sorted(m.items())) for m in masks),
                          partial, 0, 1-masks[partial][0]) is None
    assert prepared[0] in E and masks[0] == dict(enumerate(bits_of(prepared[0])))
    return {'N': N, 'scan_slots': q, 'scan_Reads': 15*q,
            'literal_operand_sum': total_length, 'literal_operand_max': maximum,
            'preparation_events': 75*N+4, 'indexed_partial_support_before_completion': 27,
            'indexed_unread_support': 1944, 'original_unread_index_1based': unread+1,
            'full_tree_vs_local_Read_checks': probes, 'partial_completion_zero_probabilities': updates,
            'completed_replay_support': 1, 'contrary_completed_bit_support': 0}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=Path(__file__).with_name('exact-results.json'))
    args = parser.parse_args()
    E0, E1 = literal_edges(0), literal_edges(1)
    geometries = []
    for i, E in enumerate((E0, E1)):
        rows, cols = index(E)
        assert len(E) == M and len(rows) == 72 and len(cols) == 243
        assert {len(v) for v in rows.values()} == {27}
        assert {len(v) for v in cols.values()} == {8}
        geometry = source.twin_geometry(i, source.literal(i))
        lc, rc, incidences, _, summary = geometry
        # All 72*243 literal incidences, including every free w3,z lift.
        assert {(a, t) for x, y in incidences for a in lc[x] for t in rc[y]} == E
        assert sum(len(lc[x])*len(rc[y]) for x, y in incidences) == M
        geometries.append(summary)
    D0, D1 = g0_witness(E0), g1_witness(E1)
    assert g0_complete(D0) == (E0, [144, 72])
    anchors_only = {(a, t) for a, t in D0 if RIGHT[t][1] == (0, 0)}
    assert len(anchors_only) == 216 and tuple(map(len, index(anchors_only))) == (72, 27)
    assert g0_complete(anchors_only) is None
    assert g1_complete(D1) == E1
    assert [len(a)*len(t) for a, t in components(D1)] == [216]*9
    # Exhaustive literal triple boundary, not a quotient stand-in.
    max_common = []
    k310 = []
    for E in (E0, E1):
        rows, _ = index(E)
        counts = Counter(len(rows[a] & rows[b] & rows[c]) for a, b, c in combinations(range(72), 3))
        max_common.append(max(counts))
        k310.append(sum(6*n*source.falling(h, 10) for h, n in counts.items()))
    assert max_common == [9, 27] and k310 == [0, 92575499208192000]
    # Every unexposed literal endpoint admits a different positive full injection.
    support_counterexamples = 0
    for E in (E0, E1):
        for side, count, ambient in ((0, 72, 128), (1, 243, 256)):
            for x in range(count):
                D = {y for y in E if y[side] != x}
                lm, rm = list(range(128)), list(range(256))
                (lm if side == 0 else rm)[x] = count
                other = transport(E, lm, rm)
                assert D <= other and other != E and len(other) == M
                assert (len({y[side] for y in other}), ambient-count) == (count, ambient-count)
                support_counterexamples += 1
    # Removing one of the two initial-stage links creates explicit alternative E.
    weakened_g0 = 0
    for t, (w, z) in enumerate(RIGHT):
        if w[0] in (0, 1) and z != (0, 0):
            second = LEFT.index(((1, 0), source.bilinear(0, (1, 0), w)))
            lost = (second, t)
            assert lost in D0
            w2 = (2, w[1], w[2])
            u = RIGHT.index((w2, (0, 1)))
            D = D0-{lost}
            other = transport(E0, range(128), swap(256, t, u))
            assert D <= other and E0 != other
            assert g0_complete(D) is None
            weakened_g0 += 1
    # Full support with 243 pairs is insufficient, even under source one alone.
    cc1 = components(E1)
    D243 = {(aa[j % 8], t) for aa, tt in cc1 for j, t in enumerate(tt)}
    assert len(D243) == 243 and tuple(map(len, index(D243))) == (72, 243)
    aa, tt = cc1[0]; bb, uu = cc1[1]
    lm, rm = swap(128, aa[0], bb[0]), list(range(256))
    for t, u in zip(tt[0::8], uu[0::8]):
        rm[t], rm[u] = rm[u], rm[t]
    other = transport(E1, lm, rm)
    assert transport(D243, lm, rm) == D243 and D243 <= other and other != E1
    # Bare double stars in G0 and G1 admit identical labelled positive data.
    stars0 = []
    for x, y in product(range(3), repeat=2):
        anchor = LEFT.index(((1, 0), (x, y)))
        tt = sorted(t for a, t in E0 if a == anchor)
        r = RIGHT.index(((x, y, -x % 3), (0, 0)))
        aa = sorted(a for a, t in E0 if t == r)
        stars0.append(([anchor]+[a for a in aa if a != anchor], [r]+[t for t in tt if t != r]))
    assert len({a for aa, tt in stars0 for a in aa}) == 72
    assert len({t for aa, tt in stars0 for t in tt}) == 243
    lm, rm = [None]*72, [None]*243
    D306 = set()
    for (aa, tt), (bb, uu) in zip(stars0, cc1):
        for a, b in zip(aa, bb): lm[a] = b
        for t, u in zip(tt, uu): rm[t] = u
        D306.update((bb[0], u) for u in uu)
        D306.update((b, uu[0]) for b in bb[1:])
    other0 = transport(E0, lm, rm)
    assert len(D306) == 306 and D306 <= other0 and D306 <= E1 and other0 != E1
    assert g1_complete(D306) is None and g1_complete(D306, both=False) == E1
    # Occupancy checks count actual 1944-edge draws. At n=q no repeat can cover.
    occupancy_results = {}
    for q in (576, 324):
        assert occupancy(q, 4) == occupancy(q, q-1) == 0
        pq, pq1 = occupancy(q, q), occupancy(q, q+1)
        assert pq == Fraction(factorial(q), M**q)
        # Independent enumeration by the one repetition/nonwitness alternative.
        independent = Fraction(factorial(q)*(q+1)*(2*M-q), 2*M**(q+1))
        assert pq1 == independent > pq > 0
        occupancy_results[str(q)] = {'n4': '0', 'nq_minus1': '0',
            'nq': str(pq), 'nq_plus1': str(pq1),
            'N_equal_q_unread_after_full_scan': 0, 'N_equal_q_plus1_reserved_unread': 1}
    # Native decoder executes the entire literal scan for each actual witness.
    native = [native_checks(E0, D0), native_checks(E1, D1)]
    result = {'sources': geometries, 'edges_each': M, 'literal_incidence_checks_each': 72*243,
              'witnesses': [{'source': 0, 'pairs': 576, 'supports': [72, 243],
                             'completion_round_sizes': [144, 72], 'recovered_edges': M},
                            {'source': 1, 'pairs': 324, 'supports': [72, 243],
                             'evidence_component_edges': sorted(sum((a,t) in D1 for a in aa for t in tt) for aa,tt in cc1),
                             'recovered_edges': M}],
              'all_literal_left_triples_each': 59640, 'maximum_common_neighbors': max_common,
              'ordered_K310_embeddings': k310,
              'unseen_endpoint_full_injection_counterexamples': support_counterexamples,
              'weakened_G0_two_link_counterexamples': weakened_g0,
              'G0_anchor_only_216_pair_support_rejection': True,
              'full_support_243_pair_noncertificate': True,
              'two_source_306_pair_double_star_noncertificate': True,
              'occupancy': occupancy_results, 'native_indexed_consumer': native,
              'all_repeat_failure': {'pairs': 1, 'supports': [1, 1], 'recognizer': 'reject',
                                     'probability_n_ge_1': '1944^(1-n)'},
              'suppliers': {'revision': BASE, 'path': 'rro-structural-recovery/exact_checks.py',
                            'sha256': hashlib.sha256(blob).hexdigest()},
              'scope': 'Exact finite validation, not a full-injection enumeration, Lean proof, optimality or originality claim.'}
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print(json.dumps({'result': 'passed', 'witness_pairs': [576,324], 'recovered_edges': [1944,1944],
                      'support_counterexamples': support_counterexamples, 'weakened_G0_counterexamples': weakened_g0,
                      'native_Reads': [x['scan_Reads'] for x in native], 'output': str(args.output)}, sort_keys=True))


if __name__ == '__main__':
    main()
