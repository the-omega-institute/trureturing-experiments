from collections import Counter
from functools import lru_cache
from itertools import combinations
from math import comb

N = 14
ALL = (1 << N) - 1
D = sum(1 << i for i in range(0, 12, 2))

def catalan(k):
    return comb(2 * k, k) // (k + 1)

def independent(s):
    return not (s & (((s << 1) & ALL) | (s >> (N - 1))))

I = [s for s in range(1 << N) if independent(s)]
W = {s: sum((-1) ** ((t ^ s).bit_count()) * catalan(12 - t.bit_count())
            for t in I if t & s == s) for s in I}
assert all(w >= 0 for w in W.values())
W = {s: w for s, w in W.items() if w > 0}
assert (len(I), len(W), sum(W.values())) == (843, 828, 208012)

# Independent count by the unique triangle adjacent to the root side.
@lru_cache(None)
def triang(a, b):
    if b - a <= 1:
        return (0,)
    out = []
    for k in range(a + 1, b):
        ears = (1 << k) if b - a == 2 else 0
        if a == 0 and b == 13 and k == 1:
            ears |= 1 << 0
        if a == 0 and b == 13 and k == 12:
            ears |= 1 << 13
        for left in triang(a, k):
            for right in triang(k, b):
                out.append(left | right | ears)
    return tuple(out)
assert Counter(triang(0, 13)) == W

def pairings(xs):
    if not xs:
        yield ()
        return
    for y in xs[1:]:
        for rest in pairings(tuple(x for x in xs[1:] if x != y)):
            yield ((xs[0], y),) + rest

# Exact recurrence for maximal cliques after adjoining one vertex.
def maximal_cliques_incremental(vs):
    cliques = {frozenset()}
    for j, (s, _, _) in enumerate(vs):
        candidates = set(cliques)
        for c in cliques:
            candidates.add(frozenset(k for k in c if vs[k][0] & s) | {j})
        cliques = {c for c in candidates if not any(c < d for d in candidates)}
    return cliques

# A separate pivoted clique enumeration, with a different recurrence.
def maximal_cliques_pivot(vs):
    nv = len(vs)
    adj = [sum(1 << j for j in range(nv)
               if i != j and vs[i][0] & vs[j][0]) for i in range(nv)]
    answer = set()
    def visit(r, p, x):
        if p == 0 and x == 0:
            answer.add(frozenset(i for i in range(nv) if r >> i & 1))
            return
        u = max((i for i in range(nv) if (p | x) >> i & 1),
                key=lambda i: (p & adj[i]).bit_count())
        candidates = p & ~adj[u]
        while candidates:
            b = candidates & -candidates
            i = b.bit_length() - 1
            visit(r | b, p & adj[i], x & adj[i])
            p ^= b
            x |= b
            candidates ^= b
    visit(0, (1 << nv) - 1, 0)
    return answer

EXPECTED = [
    (45, 27, 4290), (11, 1, 3690), (19, 6, 3690),
    (11, 1, 3690), (6, 1, 3690), (10, 2, 3690),
    (19, 6, 3690), (6, 1, 3690), (6, 1, 3690),
    (10, 2, 3690), (6, 1, 3690), (10, 2, 3690),
    (19, 6, 3690), (8, 1, 3690), (11, 1, 3690),
]
rows = []
report_rows = []
for pairs in pairings(tuple(range(0, 12, 2))):
    ms = [sum(1 << x for x in pair) for pair in pairs]
    vs = [(s, W[s], i) for i, m in enumerate(ms)
          for s in sorted(W) if s & D == m]
    assert all(sum(w for _, w, j in vs if j == i) == 4344 for i in range(3))
    mc = maximal_cliques_incremental(vs)
    assert mc == maximal_cliques_pivot(vs)
    assert all(all(vs[a][0] & vs[b][0] for a, b in combinations(c, 2)) for c in mc)
    spanning = [c for c in mc if {vs[j][2] for j in c} == {0, 1, 2}]
    row = (len(mc), len(spanning), max(sum(vs[j][1] for j in c) for c in spanning))
    rows.append(row)
    report_rows.append({
        'pairs': [list(pair) for pair in pairs],
        'vertices_per_core': [sum(j == i for _, _, j in vs) for i in range(3)],
        'all_maximal_cliques': row[0],
        'spanning_maximal_cliques': row[1],
        'maximum_spanning_weight': row[2],
    })
assert rows == EXPECTED
assert max(row[2] for row in rows) == 4290 < 4344
print('PASS: 208012 actual triangulations; all 828 positive weights agree;')
print('both clique enumerations agree on every maximal clique in all 15 partitions;')
print('every spanning clique has weight at most 4290 < 4344.')

# Optional deterministic result serialization; the proof computation above is unchanged.
if __name__ == '__main__':
    import argparse
    import json
    from pathlib import Path
    parser = argparse.ArgumentParser(description='Reproduce the supplied fourteen-gon finite table.')
    parser.add_argument('--results', type=Path)
    args = parser.parse_args()
    if args.results is not None:
        report = {
            'n': N,
            'selected_centres': list(range(0, 12, 2)),
            'independent_masks': len(I),
            'positive_complete_ear_masks': len(W),
            'actual_triangulations': sum(W.values()),
            'two_ear_weight_range': [min(w for s, w in W.items() if s.bit_count() == 2),
                                     max(w for s, w in W.items() if s.bit_count() == 2)],
            'rows': report_rows,
            'maximum_spanning_weight': max(row[2] for row in rows),
            'exact_trace_core_weight': 4344,
            'enumerations_agree': True,
        }
        args.results.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
