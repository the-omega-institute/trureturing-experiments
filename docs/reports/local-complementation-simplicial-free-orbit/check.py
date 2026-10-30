#!/usr/bin/env python3
"""Exact LC-orbit checks; standard library only, no floating-point evidence."""
from collections import defaultdict
from itertools import combinations
import json


def graph(n, edges):
    rows = [0] * n
    for u, v in edges:
        assert 0 <= u < n and 0 <= v < n and u != v
        rows[u] |= 1 << v
        rows[v] |= 1 << u
    return tuple(rows)


def edges(g):
    return [(u, v) for u in range(len(g)) for v in range(u + 1, len(g))
            if g[u] >> v & 1]


def lc(g, v):
    rows = list(g)
    for u in range(len(g)):
        if g[v] >> u & 1:
            rows[u] ^= g[v] & ~(1 << u)
    return tuple(rows)


def lc_independent(g, v):
    es = set(edges(g))
    for u, w in combinations([u for u in range(len(g)) if g[v] >> u & 1], 2):
        if (u, w) in es:
            es.remove((u, w))
        else:
            es.add((u, w))
    return graph(len(g), es)


def nonclique_witness(g, v):
    return next(((u, w) for u, w in combinations(
        [u for u in range(len(g)) if g[v] >> u & 1], 2)
        if not (g[u] >> w & 1)), None)


def connected(g):
    seen = {0}
    todo = [0]
    while todo:
        u = todo.pop()
        for v in range(len(g)):
            if g[u] >> v & 1 and v not in seen:
                seen.add(v)
                todo.append(v)
    return len(seen) == len(g)


def signatures(g):
    return tuple((g[u].bit_count(), tuple(sorted(g[v].bit_count()
                  for v in range(len(g)) if g[u] >> v & 1)))
                 for u in range(len(g)))


def isomorphism(g, h):
    """Return p with Adj_g(u,v) = Adj_h(p[u],p[v]), by exact backtracking."""
    n = len(g)
    sg, sh = signatures(g), signatures(h)
    if sorted(sg) != sorted(sh):
        return None
    p = [-1] * n
    used = set()

    def search():
        if len(used) == n:
            return tuple(p)
        choices = []
        for u in range(n):
            if p[u] != -1:
                continue
            candidates = [v for v in range(n) if v not in used and sg[u] == sh[v]
                          and all(((g[u] >> x) & 1) == ((h[v] >> p[x]) & 1)
                                  for x in range(n) if p[x] != -1)]
            if not candidates:
                return None
            choices.append((len(candidates), u, candidates))
        _, u, candidates = min(choices)
        for v in candidates:
            p[u] = v
            used.add(v)
            answer = search()
            if answer is not None:
                return answer
            used.remove(v)
            p[u] = -1
        return None

    return search()


def orbit(seed):
    representatives = [seed]
    buckets = defaultdict(list)
    buckets[tuple(sorted(signatures(seed)))].append(0)
    transitions = []
    i = 0
    while i < len(representatives):
        g = representatives[i]
        for v in range(len(g)):
            target = lc(g, v)
            key = tuple(sorted(signatures(target)))
            match = None
            for j in buckets[key]:
                p = isomorphism(target, representatives[j])
                if p is not None:
                    match = (j, p)
                    break
            if match is None:
                j = len(representatives)
                representatives.append(target)
                buckets[key].append(j)
                match = (j, tuple(range(len(g))))
            j, p = match
            transitions.append((i, v, j, p))
        i += 1
    return representatives, transitions


def verify_certificate(seed, representatives, transitions):
    n = len(seed)
    assert representatives[0] == seed
    assert len(transitions) == len(representatives) * n
    assert {(i, v) for i, v, _, _ in transitions} == {
        (i, v) for i in range(len(representatives)) for v in range(n)}
    for g in representatives:
        assert connected(g)
        assert all(nonclique_witness(g, v) is not None for v in range(n))
    for i, v, j, p in transitions:
        g, h = representatives[i], representatives[j]
        target = lc_independent(g, v)
        assert target == lc(g, v)
        assert lc(target, v) == g
        assert sorted(p) == list(range(n))
        assert all(((target[u] >> w) & 1) == ((h[p[u]] >> p[w]) & 1)
                   for u in range(n) for w in range(n))


def wheel():
    return graph(6, [(i, (i + 1) % 5) for i in range(5)] + [(i, 5) for i in range(5)])


def main():
    # Source IV.A: local complementation maps K_4 to a star at its chosen vertex.
    k4 = graph(4, combinations(range(4), 2))
    star = graph(4, [(0, i) for i in range(1, 4)])
    assert lc(k4, 0) == star == lc_independent(k4, 0)
    assert all(nonclique_witness(star, v) is None for v in [1, 2, 3])
    print('FIDELITY_K4_TO_STAR_OK')

    # Source Fig. 1 / Appendix E class 19: the six-vertex wheel orbit.
    w = wheel()
    reps6, trans6 = orbit(w)
    assert len(reps6) == 2
    verify_certificate(w, reps6, trans6)
    source_rep = graph(6, [(0, 1), (0, 2), (0, 5), (1, 2), (1, 4),
                           (2, 3), (3, 4), (3, 5), (4, 5)])
    assert any(isomorphism(source_rep, h) is not None for h in reps6)
    print('FIDELITY_CLASS_19_ORBIT_OK representatives=2 transitions=12')

    seed = graph(12, edges(w) + [(u + 6, v + 6) for u, v in edges(w)] + [(5, 11)])
    representatives, transitions = orbit(seed)
    assert len(representatives) == 54
    assert len(transitions) == 648
    assert len(edges(seed)) == 21
    verify_certificate(seed, representatives, transitions)
    print('CANDIDATE_ORBIT_OK vertices=12 edges=21 representatives=54 transitions=648')
    print(json.dumps({'seed_edges': edges(seed),
                      'representatives': representatives,
                      'transitions': transitions,
                      'nonclique_witnesses': [
                          [nonclique_witness(g, v) for v in range(12)]
                          for g in representatives]}, separators=(',', ':')))
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
