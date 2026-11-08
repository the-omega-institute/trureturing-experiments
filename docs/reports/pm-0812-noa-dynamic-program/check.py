#!/usr/bin/env python3
"""Exact finite checks for the all-nOA closure candidate; not a proof.

The source is arXiv:0812.3072v1, Definition 7 and equation (6).
Greechie inputs are mathematical data transcribed from that chapter and
arXiv:quant-ph/0009038v3; no third-party program is used or copied.
"""
from itertools import product


class OL:
    def __init__(self, name, labels, le, comp):
        self.name, self.labels, self.le, self.comp = name, labels, le, comp
        self.m = len(labels)
        n = self.m
        self.meet = [[0] * n for _ in range(n)]
        self.join = [[0] * n for _ in range(n)]
        for x, y in product(range(n), repeat=2):
            lower = [z for z in range(n) if le[z][x] and le[z][y]]
            upper = [z for z in range(n) if le[x][z] and le[y][z]]
            greatest = [z for z in lower if all(le[w][z] for w in lower)]
            least = [z for z in upper if all(le[z][w] for w in upper)]
            assert len(greatest) == len(least) == 1, (name, x, y)
            self.meet[x][y], self.join[x][y] = greatest[0], least[0]
        self.bot = next(z for z in range(n) if all(le[z]))
        self.top = next(z for z in range(n) if all(le[x][z] for x in range(n)))
        for x in range(n):
            assert le[x][x] and comp[comp[x]] == x
            assert self.meet[x][comp[x]] == self.bot
            assert self.join[x][comp[x]] == self.top
        for x, y, z in product(range(n), repeat=3):
            assert not (le[x][y] and le[y][z]) or le[x][z]
            assert not (le[x][y] and le[y][x]) or x == y
            assert not le[x][y] or le[comp[y]][comp[x]]
        height_at = [0] * n
        for x in sorted(range(n), key=lambda x: sum(le[y][x] for y in range(n))):
            height_at[x] = max([0] + [height_at[y] + 1 for y in range(n)
                                     if y != x and le[y][x]])
        self.h = max(height_at)
        self.imp = [[self.join[comp[x]][self.meet[x][c]] for c in range(n)]
                    for x in range(n)]
        self.oml = all(not le[x][y] or
                       self.join[x][self.meet[comp[x]][y]] == y
                       for x, y in product(range(n), repeat=2))


def from_covers(name, labels, covers, comp):
    n = len(labels)
    le = [[x == y for y in range(n)] for x in range(n)]
    for x, y in covers:
        le[x][y] = True
    for z, x, y in product(range(n), repeat=3):
        le[x][y] |= le[x][z] and le[z][y]
    return OL(name, labels, le, comp)


def boolean(d):
    n = 1 << d
    return OL('Boolean-' + str(d), [str(i) for i in range(n)],
              [[x & y == x for y in range(n)] for x in range(n)],
              [x ^ (n - 1) for x in range(n)])


def greechie(name, text):
    blocks = [tuple(b) for b in text.rstrip('.').split(',')]
    assert all(len(b) == 3 and len(set(b)) == 3 for b in blocks)
    atoms = sorted(set(''.join(''.join(b) for b in blocks)))
    k = len(atoms)
    ai = {a: i + 1 for i, a in enumerate(atoms)}
    n = 2 * k + 2
    top = n - 1
    labels = ['0'] + atoms + [a + "'" for a in atoms] + ['1']
    comp = [top] + [i + k for i in range(1, k + 1)] + list(range(1, k + 1)) + [0]
    covers = [(0, i) for i in range(1, k + 1)]
    covers += [(i + k, top) for i in range(1, k + 1)]
    for b in blocks:
        covers += [(ai[x], ai[y] + k) for x in b for y in b if x != y]
    return from_covers(name, labels, covers, comp)


# Algorithm under test: simultaneous whole-matrix pivot updates.
def initial(L, c):
    A, C, J, I = L.imp, L.comp, L.join, L.meet
    return tuple(tuple(J[I[A[x][c]][A[y][c]]][I[A[C[x]][c]][A[C[y]][c]]]
                       for y in range(L.m)) for x in range(L.m))


def pivot(L, M, v):
    I, J = L.meet, L.join
    return tuple(tuple(J[M[x][y]][I[M[x][v]][M[y][v]]]
                       for y in range(L.m)) for x in range(L.m))


def closure(L, c):
    M, word, sweeps, visits = initial(L, c), [], 0, 0
    N = L.h * L.m ** 2
    while True:
        before = M
        sweeps += 1
        assert sweeps <= N + 1
        for v in range(L.m):
            visits += 1
            Q = pivot(L, M, v)
            assert all(L.le[M[x][y]][Q[x][y]]
                       for x, y in product(range(L.m), repeat=2))
            if Q != M:
                word.append(v)
                assert len(word) <= N
                M = Q
        if M == before:
            break
    assert visits <= L.m * (N + 1)
    assert all(pivot(L, M, v) == M for v in range(L.m))
    replay = initial(L, c)
    for v in word:
        replay = pivot(L, replay, v)
    assert replay == M
    return M, tuple(word), sweeps, visits


def criterion(L, c, M):
    return all(L.le[L.meet[L.imp[x][c]][M[x][y]]][L.imp[y][c]]
               for x, y in product(range(L.m), repeat=2))


# Independent scalar evaluator: literal source recursion, not calling pivot,
# initial, or the algorithm's precomputed implication table.
def source_term(L, c, tail, x, y):
    I, J, C = L.meet, L.join, L.comp
    if not tail:
        xc = J[C[x]][I[x][c]]
        yc = J[C[y]][I[y][c]]
        xpc = J[C[C[x]]][I[C[x]][c]]
        ypc = J[C[C[y]]][I[C[y]][c]]
        return J[I[xc][yc]][I[xpc][ypc]]
    v, rest = tail[-1], tail[:-1]
    return J[source_term(L, c, rest, x, y)][
        I[source_term(L, c, rest, x, v)][source_term(L, c, rest, y, v)]]


def brute_law(L, t, terminals):
    valuations = 0
    # Each negative result returns an actual assignment, not a failed bound.
    for c in range(L.m):
        previous_prefix = None
        prefix_values = None
        for tail in product(range(L.m), repeat=t):
            if t == 2 and tail[:-1] != previous_prefix:
                previous_prefix = tail[:-1]
                # Partial evaluation of the literal scalar recurrence: all
                # three uses below have this SAME fixed prefix assignment.
                prefix_values = tuple(tuple(source_term(L, c, tail[:-1], x, y)
                                            for y in range(L.m)) for x in range(L.m))
            M = initial(L, c)
            for v in tail:
                M = pivot(L, M, v)
            for x, y in product(range(L.m), repeat=2):
                valuations += 1
                if t == 2:
                    v = tail[-1]
                    s = L.join[prefix_values[x][y]][
                        L.meet[prefix_values[x][v]][prefix_values[y][v]]]
                else:
                    s = source_term(L, c, tail, x, y)
                assert s == M[x][y], (L.name, c, tail, x, y, s, M[x][y])
                assert L.le[s][terminals[c][0][x][y]]
                lhs = L.meet[L.imp[x][c]][s]
                if not L.le[lhs][L.imp[y][c]]:
                    assert not criterion(L, c, M)
                    return False, valuations, (x, y, c, tail)
    return True, valuations, None


def main():
    o6 = from_covers('O6', ['0', 'a', 'b', "a'", "b'", '1'],
                     [(0, 1), (1, 2), (2, 5), (0, 4), (4, 3), (3, 5)],
                     [5, 3, 4, 1, 2, 0])
    mo2 = from_covers('MO2', ['0', 'a', "a'", 'b', "b'", '1'],
                      [(0, i) for i in range(1, 5)] + [(i, 5) for i in range(1, 5)],
                      [5, 2, 1, 4, 3, 0])
    try:
        greechie('chapter-caption-literal', '123,345,567,789,9AB,BC1,BD5')
    except AssertionError as e:
        print('SOURCE_DATA_BOUNDARY: literal chapter caption is not a lattice:', e, flush=True)
    else:
        raise AssertionError('Expected the chapter caption lattice defect')
    fixtures = [
        (boolean(0), True, [True, True, True]),
        (boolean(1), True, [True, True, True]),
        (boolean(2), True, [True, True, True]),
        (boolean(3), True, [True, True, True]),
        (o6, False, [False, False, False]),
        (mo2, True, [True, True, True]),
        (greechie('L28-13-7-OMLp-oa3f', '123,345,567,789,9AB,BC1,6DC'),
         False, [False, False, False]),
        (greechie('L36-17-10-oa3p4f', '8AH,89G,79F,6AE,59C,4AD,38B,245,167,123'),
         False, [True, False, False]),
        (greechie('L46-7-22-13-oa4p5f', 'DEM,BEI,AEJ,9CL,8CK,7CF,6DH,5DG,49B,38A,246,135,127'),
         False, [True, True, False]),
        (greechie('L46-9-22-13-oa4p5f', 'CDM,BFG,AEH,9FK,8EL,7EI,6FJ,5BD,4AC,257,239,146,138'),
         False, [True, True, False]),
    ]
    for L, expected_all, expected_small in fixtures:
        terminals = [closure(L, c) for c in range(L.m)]
        all_pass = all(criterion(L, c, r[0]) for c, r in enumerate(terminals))
        assert all_pass == expected_all, L.name
        print('LATTICE', L.name, 'm=' + str(L.m), 'h=' + str(L.h),
              'orthomodular=' + str(L.oml), 'all_nOA=' + str(all_pass), flush=True)
        print('CLOSURE', 'max_changes=' + str(max(len(r[1]) for r in terminals)),
              'max_sweeps=' + str(max(r[2] for r in terminals)),
              'word_bound=' + str(L.h * L.m ** 2), flush=True)
        for t, expected in enumerate(expected_small):
            actual, count, w = brute_law(L, t, terminals)
            assert actual == expected, (L.name, t, actual)
            print('SOURCE', str(t + 3) + 'OA', 'holds=' + str(actual),
                  'valuations=' + str(count),
                  'witness=' + str(None if w is None else
                                  ([L.labels[w[i]] for i in range(3)],
                                   [L.labels[v] for v in w[3]])), flush=True)
        if not all_pass:
            c = next(c for c, r in enumerate(terminals) if not criterion(L, c, r[0]))
            R, w, _, _ = terminals[c]
            x, y = next((x, y) for x, y in product(range(L.m), repeat=2)
                        if not L.le[L.meet[L.imp[x][c]][R[x][y]]][L.imp[y][c]])
            value = source_term(L, c, w, x, y)
            assert value == R[x][y]
            assert not L.le[L.meet[L.imp[x][c]][value]][L.imp[y][c]]
            # Padding by a terminal no-op keeps the same real source assignment.
            assert source_term(L, c, w + (L.bot,), x, y) == value
            print('CERTIFICATE', 'n=' + str(3 + len(w)),
                  'x,y,c=' + str([L.labels[z] for z in (x, y, c)]),
                  'tail=' + str([L.labels[v] for v in w]), flush=True)
    print('FIDELITY_OK: O6 fails 3OA; source 13-7 fails 3OA; L36 passes 3OA and fails 4OA; both L46 examples pass 4OA and fail 5OA')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
