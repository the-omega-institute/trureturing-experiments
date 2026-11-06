#!/usr/bin/env python3
"""Exact incidence, completed-placement and GF(2) checks for arXiv:2609.29149v1.

All counts, data and linear algebra use integers. Finite checks do not prove the
universal characterization. No third-party program source is used.
"""
from itertools import combinations, product, permutations
from collections import Counter


def subsets(xs, nonempty=False):
    xs = tuple(xs)
    for r in range(int(nonempty), len(xs) + 1):
        yield from (frozenset(s) for s in combinations(xs, r))


def design(A, k, lam):
    v = len(A)
    assert all(i not in A[i] and len(A[i]) == k for i in range(v))
    assert all(sum(t in row for row in A) == k for t in range(v))
    assert all(len(A[i] & A[j]) == lam for i in range(v) for j in range(i))


def blocked(A, H, F, i):
    return frozenset(t for t in A[i] if t in F or sum(t in A[j] for j in H) >= 2)


def core_value(A, H, F):
    return min(len(blocked(A, H, F, i)) for i in H)


def peel(A, F):
    H, result = set(F), 0
    while H:
        i = min(H, key=lambda i: (len(blocked(A, H, F, i)), i))
        result = max(result, len(blocked(A, H, F, i)))
        H.remove(i)
    return result


def gamma(A, p, check_peeling=True):
    answer, arg = -1, None
    for F in subsets(range(len(A)), True):
        if len(F) > p:
            continue
        values = [(core_value(A, H, F), H) for H in subsets(F, True)]
        value, H = max(values, key=lambda x: x[0])
        if check_peeling:
            assert peel(A, F) == value
        if value > answer:
            answer, arg = value, (F, H)
    assert answer >= 0
    return answer, arg


def gamma_bounded(A, p):
    # This routine avoids enumerating irrelevant subsets when v = 40.
    answer, arg, count = -1, None, 0
    for r in range(1, p + 1):
        for ft in combinations(range(len(A)), r):
            F = frozenset(ft)
            best, best_H = -1, None
            for H in subsets(F, True):
                x = core_value(A, H, F)
                if x > best:
                    best, best_H = x, H
            assert peel(A, F) == best
            count += 1
            if best > answer:
                answer, arg = best, (F, best_H)
    return answer, arg, count


def decompose(counts, N, q):
    """Degree sequence in [0,N], sum qN -> N distinct-target q-subsets."""
    counts = dict(counts)
    assert sum(counts.values()) == q * N
    assert all(0 <= a <= N for a in counts.values())
    result = []
    for remaining in range(N, 0, -1):
        selected = sorted(counts, key=lambda t: (-counts[t], t))[:q]
        assert len(selected) == q and all(counts[t] > 0 for t in selected)
        selected = frozenset(selected)
        result.append(selected)
        for t in selected:
            counts[t] -= 1
        assert all(0 <= a <= remaining - 1 for a in counts.values())
        assert sum(counts.values()) == q * (remaining - 1)
    assert not any(counts.values())
    return result


def complete_witness(A, q, H, F):
    """Build the claimed exactly 2k² objects per primary and all complete parities."""
    v, k = len(A), len(A[0])
    assert H and H <= F and all(core_value(A, H, F) >= q for _ in [0])
    M, R = 2 * k * k, 2 * q * k
    primary, selected, marked = [], [], set()
    def add(i, T, mark=False):
        x = len(primary)
        primary.append(i)
        selected.append(frozenset(T))
        if mark:
            marked.add(x)
    if q == 0:
        for i in range(v):
            for a in range(M):
                add(i, (), i == min(H) and a == 0)
    else:
        for i in range(v):
            B = sorted(blocked(A, H, F, i)) if i in H else []
            c = len(B)
            if i in H:
                assert q <= c <= k
                for s in range(c):
                    T = {B[(s + r) % c] for r in range(q)}
                    assert len(T) == q
                    add(i, T, True)
                    add(i, T, True)
            N = M - 2 * c
            counts = {t: R - (2 * q if t in B else 0) for t in A[i]}
            for T in decompose(counts, N, q):
                add(i, T)
    assert len(primary) == v * M
    edge = {(i, t): [] for i in range(v) for t in A[i]}
    for x, T in enumerate(selected):
        for t in T:
            edge[primary[x], t].append(x)
    groups = []
    for t in range(v):
        src = [i for i in range(v) if t in A[i]]
        assert len(src) == k
        assert all(len(edge[i, t]) == R for i in src)
        slots = [{} for _ in range(R)]
        special = {i: [x for x in edge[i, t] if x in marked] for i in src}
        if t not in F:
            active = [i for i in src if special[i]]
            if active:
                assert len(active) >= 2
                assert set(active) == set(src) & set(H)
                assert all(len(special[i]) == 2 * q for i in active)
                s = 0
                # A 2-cycle uses two parallel pairs of distinct replica identities.
                for _ in range(q):
                    for a in range(len(active)):
                        i, j = active[a], active[(a + 1) % len(active)]
                        slots[s][i], slots[s][j] = special[i].pop(), special[j].pop()
                        s += 1
                assert s == q * len(active) <= R
        else:
            for i in src:
                for s, x in enumerate(special[i]):
                    slots[s][i] = x
                special[i].clear()
        assert not any(special.values())
        for i in src:
            unmarked = [x for x in edge[i, t] if x not in marked]
            free = [s for s in range(R) if i not in slots[s]]
            assert len(free) == len(unmarked)
            for s, x in zip(free, unmarked):
                slots[s][i] = x
        groups.append([tuple(slot[i] for i in src) for slot in slots])
    # Independent legality and two-assignment checks on the constructed objects.
    assert Counter(primary) == {i: M for i in range(v)}
    for x, T in enumerate(selected):
        assert len(T) == q and T <= A[primary[x]]
    for t in range(v):
        received = Counter(x for G in groups[t] for x in G)
        assert received == Counter({x: 1 for x, T in enumerate(selected) if t in T})
        for G in groups[t]:
            assert len(G) == k and len(set(G)) == k
            assert Counter(primary[x] for x in G) == {i: 1 for i in range(v) if t in A[i]}
            if t not in F:
                assert sum(x in marked for x in G) in (0, 2)
    assert marked and all(primary[x] in F for x in marked)
    # data1 = 0; data2 = the nonempty marked support. Surviving observations agree.
    assert all(x not in marked for x in range(len(primary)) if primary[x] not in F)
    return primary, selected, groups, marked


def rank2(rows):
    basis = {}
    for row in rows:
        while row:
            i = row.bit_length() - 1
            if i in basis:
                row ^= basis[i]
            else:
                basis[i] = row
                break
    return len(basis)


def lost_rank(primary, groups, F):
    lost = [x for x, i in enumerate(primary) if i in F]
    idx = {x: j for j, x in enumerate(lost)}
    rows = [sum(1 << idx[x] for x in G if x in idx)
            for t, gs in enumerate(groups) if t not in F for G in gs]
    return rank2(rows), len(lost)


def regular_matrices(v):
    for k in range(1, v):
        options = [list(map(frozenset, combinations([t for t in range(v) if t != i], k)))
                   for i in range(v)]
        for A in product(*options):
            if all(sum(t in row for row in A) == k for t in range(v)):
                yield A


def exhaustive_regular_checks():
    matrices = cores = witnesses = 0
    for v in range(2, 5):
        for A in regular_matrices(v):
            matrices += 1
            for F in subsets(range(v), True):
                best = max(core_value(A, H, F) for H in subsets(F, True))
                assert best == peel(A, F)
                for H in subsets(F, True):
                    cores += 1
                    for q in range(core_value(A, H, F) + 1):
                        complete_witness(A, q, H, F)
                        witnesses += 1
    print(f'REGULAR_EXHAUSTIVE v=2..4 matrices={matrices} nonempty_cores={cores} witnesses={witnesses}')


def brute_small_placements():
    # Independently enumerate ALL choices and legal parity partitions with two objects
    # per primary, for all zero-diagonal regular representatives on 2 or 3 nodes.
    checked = 0
    for v in (2, 3):
        for A in regular_matrices(v):
            k = len(A[0])
            primary = [i for i in range(v) for _ in range(2)]
            for q in range(k + 1):
                options = [tuple(map(frozenset, combinations(A[i], q))) for i in primary]
                bad = {F: False for F in subsets(range(v), True)}
                for selected in product(*options):
                    possibilities = []
                    for t in range(v):
                        lists = [[x for x, i0 in enumerate(primary) if i0 == i and t in selected[x]]
                                 for i in range(v) if t in A[i]]
                        R = len(lists[0])
                        if any(len(xs) != R for xs in lists):
                            break
                        # Anchor the first primary's order: quotient by parity-slot permutation.
                        local = []
                        for rest in product(*(tuple(permutations(xs)) for xs in lists[1:])):
                            local.append(list(zip(lists[0], *rest)))
                        possibilities.append(local)
                    else:
                        for groups in product(*possibilities):
                            checked += 1
                            for F in bad:
                                r, n = lost_rank(primary, groups, F)
                                bad[F] |= r < n
                for F, b in bad.items():
                    obstructed = any(core_value(A, H, F) >= q for H in subsets(F, True))
                    assert b == obstructed, (A, q, F, b, obstructed)
    print(f'BRUTE_LEGAL_PLACEMENTS v=2..3 objects_per_primary=2 checked={checked} identity=OK')


def source_fidelity():
    # Exact matrix in Appendix B: complement of the Fano plane, SBIBD(7,4,2).
    rows = ('0111001', '0010111', '0101110', '1110100', '1100011', '1001101', '1011010')
    A = tuple(frozenset(t for t, c in enumerate(row) if c == '1') for row in rows)
    design(A, 4, 2)
    assert A[0] == {1, 2, 3, 6} and A[3] == {0, 1, 2, 4}
    assert A[0] & A[3] == {1, 2}
    assert gamma_bounded(A, 1)[0] == 0
    g, (F, H), count = gamma_bounded(A, 2)
    assert g == 3
    # A completed ambiguity at q=3 and exhaustive primary-type peeling at q=4.
    P = complete_witness(A, 3, frozenset({0, 3}), frozenset({0, 3}))
    r, n = lost_rank(P[0], P[2], {0, 3})
    assert r < n
    assert all(peel(A, frozenset(F)) < 4 for s in (1, 2) for F in combinations(range(7), s))
    print(f'SOURCE_FIDELITY Appendix_B SBIBD(7,4,2) p=2 Gamma={g} threshold=4 failure_sets={count} q3_rank={r}/{n} objects={len(P[0])}')
    return A


def projective_three_space(rho):
    assert rho == 3  # Prime-field arithmetic; no unsupported prime-power implementation.
    pts = [x for x in product(range(rho), repeat=4) if any(x) and next(a for a in x if a) == 1]
    special = [(1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0)]
    pts = special + [x for x in pts if x not in special]
    dot = lambda x, y: sum(a * b for a, b in zip(x, y)) % rho
    planes = list(pts)
    chosen = [(1, 1, 0, 0), (0, 1, 1, 0), (1, 0, 1, 0)]
    match = {planes.index(y): i for i, y in enumerate(chosen)}
    def augment(i, seen):
        for j, y in enumerate(planes):
            if j in seen or not dot(pts[i], y):
                continue
            seen.add(j)
            if j not in match or (match[j] >= 3 and augment(match[j], seen)):
                match[j] = i
                return True
        return False
    for i in range(3, len(pts)):
        assert augment(i, set())
    host_plane = {i: planes[j] for j, i in match.items()}
    A = tuple(frozenset(t for t in range(len(pts)) if not dot(x, host_plane[t])) for x in pts)
    design(A, 13, 4)
    assert [[int(t in A[i]) for t in range(3)] for i in range(3)] == [[0, 1, 0], [0, 0, 1], [1, 0, 0]]
    assert all(any(sum(t in A[i] for i in F) == 3 for t in range(len(A))) for F in combinations(range(len(A)), 3))
    g, (F, H), count = gamma_bounded(A, 3)
    assert g == 8 and g + 1 == 2 * rho + 3
    P = complete_witness(A, g, H, F)
    print(f'PG3_EXACT rho=3 v=40 k=13 lambda=4 p=3 no_separated_triples=OK failure_sets={count} Gamma={g} predicted_threshold={g+1} objects={len(P[0])} objects_per_primary=338 alive_XOR_ambiguity=OK')


def main():
    # Degree completion lemma on all small integer count vectors, independent of designs.
    count = 0
    for k in range(1, 6):
        for N in range(5):
            for q in range(k + 1):
                for a in product(range(N + 1), repeat=k):
                    if sum(a) == q * N:
                        ds = decompose(dict(enumerate(a)), N, q)
                        assert Counter(t for T in ds for t in T) == Counter({t: x for t, x in enumerate(a) if x})
                        count += 1
    print(f'DEGREE_COMPLETION k=1..5 N=0..4 vectors={count} OK')
    exhaustive_regular_checks()
    brute_small_placements()
    source_fidelity()
    projective_three_space(3)
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
