#!/usr/bin/env python3
"""Exact literal-face checks for Rather, arXiv:2605.30781v1, Conjecture 3.7.

Original implementation. A face is an integer bit mask on a fixed ambient set.
No conjectured recurrence is used to compute an orbit. Finite checks are not a
proof of the universal conjecture. The known P4 exception is reported separately.
"""
from itertools import combinations, permutations


def submask(a, b):
    return a & b == a


def facets(faces, vertex_count):
    """For a lower set, maximality is equivalent to no one-vertex extension."""
    full = (1 << vertex_count) - 1
    return frozenset(s for s in faces if not any(
        (s | (1 << v)) in faces
        for v in range(vertex_count) if not s & (1 << v)))


def lower_closure(generators, vertex_count):
    return frozenset(s for s in range(1 << vertex_count)
                     if any(submask(s, f) for f in generators))


def nf(faces, vertex_count):
    """Literal squarefree-monomial membership: avoid every facet generator."""
    fs = facets(faces, vertex_count)
    return frozenset(s for s in range(1 << vertex_count)
                     if all(not submask(f, s) for f in fs))


def graph_edges(n, m):
    return frozenset([(1 << a) | (1 << b)
                      for block in [range(n), range(n, n + m)]
                      for a, b in combinations(block, 2)] + [1 | (1 << n)])


def c_state(n, m, r):
    a, b, e = (1 << n) - 1, ((1 << m) - 1) << n, 1 | (1 << n)

    def in_u(s):
        return not any(submask(f, s) for f in [a, b, e])

    return frozenset(s for s in range(1 << (n + m)) if
                     s.bit_count() <= r or
                     (in_u(s) and s.bit_count() <= r + 1) or
                     (submask(s, a) and s.bit_count() <= r + 2 <= n) or
                     (submask(s, b) and s.bit_count() <= r + 2 <= m))


def predicted_c_facets(n, m, r):
    a, b, e = (1 << n) - 1, ((1 << m) - 1) << n, 1 | (1 << n)
    result = set()
    for s in range(1 << (n + m)):
        size = s.bit_count()
        in_u = not any(submask(f, s) for f in [a, b, e])
        pure = submask(s, a) or submask(s, b)
        if ((size == r and not in_u) or
                (size == r + 1 and in_u and not pure) or
                (size == r + 2 and pure)):
            result.add(s)
    return frozenset(result)


def has_triangle(faces, vertex_count):
    return any(all(((1 << a) | (1 << b)) in faces
                   for a, b in combinations(t, 2))
               for t in combinations(range(vertex_count), 3))


def image(s, perm):
    return sum(1 << perm[i] for i in range(len(perm)) if s & (1 << i))


def main():
    recurrence_checks = 0
    facet_checks = 0
    tested = 0
    # All unordered pairs in this complete finite square, including n=2 and n!=m.
    # Reversing the blocks is an explicit relabelling, so each ordered pair is covered.
    for n in range(2, 7):
        for m in range(n, 7):
            N = n + m
            a, b, e, full = (1 << n) - 1, ((1 << m) - 1) << n, 1 | (1 << n), (1 << N) - 1
            g = lower_closure(graph_edges(n, m), N)
            orbit = [g]
            for _ in range(N + 2):
                orbit.append(nf(orbit[-1], N))
            # Fidelity: Theorems 3.3 and 3.5, independently computed from all faces.
            cross = frozenset((1 << i) | (1 << j)
                              for i in range(n) for j in range(n, N)
                              if (1 << i) | (1 << j) != e)
            assert facets(orbit[1], N) == cross, (n, m, 'source Thm 3.3')
            assert facets(orbit[2], N) == frozenset([a, b, e]), (n, m, 'source Thm 3.5')
            u = frozenset(s for s in range(1 << N)
                          if not any(submask(f, s) for f in [a, b, e]))
            assert orbit[3] == u
            uf = frozenset(full ^ ((1 << i) | (1 << j))
                           for i in range(n) for j in range(n, N) if i == 0 or j == n)
            assert facets(u, N) == uf
            if n == m == 2:
                witnesses = [p for p in permutations(range(N))
                             if frozenset(image(s, p) for s in g) == orbit[1]]
                assert witnesses and orbit[1] != g and orbit[2] == g
                print(f'KNOWN_EXCEPTION n=2 m=2 first_isomorphic_return=1 labeled_return=2 permutation={witnesses[0]}')
                continue
            # Check the proposed full orbit, its entry, every recurrence, and terminal step.
            assert orbit[5] == c_state(n, m, N - 3), (n, m, 'D^5 base')
            for r in range(1, N - 2):
                cr = c_state(n, m, r)
                assert facets(cr, N) == predicted_c_facets(n, m, r), (n, m, r, 'facets')
                facet_checks += 1
                target = g if r == 1 else c_state(n, m, r - 1)
                assert nf(cr, N) == target, (n, m, r, 'recurrence/terminal')
                recurrence_checks += 1
                assert orbit[N + 2 - r] == cr, (n, m, r, 'actual orbit')
            assert orbit[N + 2] == g
            # Isomorphism invariants exclude ALL earlier returns, not just labelled ones.
            assert has_triangle(g, N) and not has_triangle(orbit[1], N)
            assert all(any(s.bit_count() >= 3 for s in orbit[k]) for k in range(2, N + 2))
            assert max(s.bit_count() for s in g) == 2
            tested += 1
            print(f'PASS n={n} m={m} first_isomorphic_return={N + 2} states_checked={N + 3} fidelity=Thm3.3,Thm3.5')
    print(f'SUMMARY nonexceptional_pairs={tested} facet_classifications={facet_checks} recurrence_and_terminal_checks={recurrence_checks} literal_vertex_limit=12')
    print('SCOPE finite exact checks only; no universal or Lean proof claimed')
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
