#!/usr/bin/env python3
"""Exact source-fidelity and graph/shear/layout checks. Python 3, standard library only."""
from collections import deque
from fractions import Fraction
from itertools import combinations, permutations, product


def transpose(a):
    return [list(row) for row in zip(*a)]


def matmul(a, b):
    bt = transpose(b)
    return [[sum(x * y for x, y in zip(row, col)) for col in bt] for row in a]


def omega(n):
    return [[(1 if i % 2 == 0 else -1) if i // 2 == j // 2 and i != j else 0
             for j in range(2 * n)] for i in range(2 * n)]


def shear(adj):
    n = len(adj)
    s = [[int(i == j) for j in range(2 * n)] for i in range(2 * n)]
    for v in range(n):
        for w in range(n):
            s[2 * v + 1][2 * w] = adj[v][w]
    return s


def supports(s):
    n = len(s) // 2
    return [sum(1 << a for a in range(n)
                if any(s[2 * o + r][2 * a + c] != 0 for r in range(2) for c in range(2)))
            for o in range(n)]


def source_profile(supp, order):
    loaded = 0
    costs = []
    for k, o in enumerate(order, 1):
        loaded |= supp[o]
        costs.append(loaded.bit_count() - k + 1)
    return costs


def rank(a):
    a = [[Fraction(x) for x in row] for row in a]
    r = 0
    for c in range(len(a[0]) if a else 0):
        pivot = next((i for i in range(r, len(a)) if a[i][c]), None)
        if pivot is None:
            continue
        a[r], a[pivot] = a[pivot], a[r]
        q = a[r][c]
        a[r] = [x / q for x in a[r]]
        for i in range(r + 1, len(a)):
            q = a[i][c]
            if q:
                a[i] = [x - q * y for x, y in zip(a[i], a[r])]
        r += 1
    return r


def fidelity_rank(s, order):
    supp = supports(s)
    loaded = 0
    for k, o in enumerate(order, 1):
        loaded |= supp[o]
        rows = [2 * v + r for v in order[k - 1:] for r in range(2)]
        cols = [2 * v + r for v in range(len(order)) if loaded >> v & 1 for r in range(2)]
        observed = rank([[s[i][j] for j in cols] for i in rows])
        assert observed == 2 * (loaded.bit_count() - k + 1), (s, order, k, observed)
    return len(order)


def neighbours(adj):
    return [sum(1 << j for j, a in enumerate(row) if a) for row in adj]


def boundary(nb, prefix):
    outside = ((1 << len(nb)) - 1) ^ prefix
    return sum(1 << v for v in range(len(nb)) if prefix >> v & 1 and nb[v] & outside)


def vertex_separation(nb, order):
    prefix = 0
    result = 0
    for v in order:
        prefix |= 1 << v
        result = max(result, boundary(nb, prefix).bit_count())
    return result


def valid_decomposition(nb, bags):
    n = len(nb)
    assert set.union(*bags) == set(range(n))
    for v in range(n):
        positions = [i for i, b in enumerate(bags) if v in b]
        assert positions == list(range(min(positions), max(positions) + 1))
        for w in range(v + 1, n):
            if nb[v] >> w & 1:
                assert any(v in b and w in b for b in bags)


def order_bags(nb, order):
    prefix = 0
    bags = []
    for v in order:
        b = boundary(nb, prefix) | (1 << v)
        bags.append({w for w in range(len(nb)) if b >> w & 1})
        prefix |= 1 << v
    return bags


def exact_pathwidth(nb):
    """Independent interval/bag search, not a vertex-separation recurrence.

    Each vertex is unintroduced, live, or forgotten. Introduction and forgetting
    produce its interval of bags. Forgetting requires every neighbour already
    introduced, ensuring edge coverage. Search all such paths under a bag cap.
    Any path decomposition can be refined to these introduce/forget moves.
    """
    n = len(nb)
    full = (1 << n) - 1
    for cap in range(1, n + 1):
        seen = {(0, 0)}
        todo = deque([(0, 0)])
        while todo:
            introduced, live = todo.popleft()
            if introduced == full and live == 0:
                return cap - 1
            for v in range(n):
                bit = 1 << v
                if not introduced & bit and live.bit_count() < cap:
                    state = (introduced | bit, live | bit)
                elif live & bit and not (nb[v] & (full ^ introduced)):
                    state = (introduced, live ^ bit)
                else:
                    continue
                if state not in seen:
                    seen.add(state)
                    todo.append(state)
        # Increasing cap preserves all previously feasible bag sequences.
    raise AssertionError('No decomposition found')


def adjacency(n, edge_mask):
    a = [[0] * n for _ in range(n)]
    for i, (v, w) in enumerate(combinations(range(n), 2)):
        if edge_mask >> i & 1:
            a[v][w] = a[w][v] = 1
    return a


def main():
    graphs = orders = rank_cuts = 0
    per_size = []
    for n in range(1, 6):
        size_orders = 0
        for mask in range(1 << (n * (n - 1) // 2)):
            adj = adjacency(n, mask)
            nb = neighbours(adj)
            s = shear(adj)
            assert matmul(matmul(transpose(s), omega(n)), s) == omega(n)
            supp = supports(s)
            assert supp == [nb[v] | (1 << v) for v in range(n)]
            best_cost = n + 1
            best_vs = n + 1
            for order in permutations(range(n)):
                cost = max(source_profile(supp, order))
                vs = vertex_separation(nb, order[::-1])
                assert cost == 1 + vs, (n, mask, order, cost, vs)
                bags = order_bags(nb, order[::-1])
                valid_decomposition(nb, bags)
                assert max(map(len, bags)) - 1 == vs
                best_cost = min(best_cost, cost)
                best_vs = min(best_vs, vs)
                if n <= 4:
                    rank_cuts += fidelity_rank(s, order)
                size_orders += 1
            pw = exact_pathwidth(nb)
            assert best_vs == pw and best_cost == 1 + pw, (n, mask, best_cost, pw)
            if not mask:
                assert best_cost == 1
            if mask == (1 << (n * (n - 1) // 2)) - 1:
                assert best_cost == n
            graphs += 1
        orders += size_orders
        per_size.append((n, 1 << (n * (n - 1) // 2), size_orders))
    print('All labelled simple graphs: (vertices, graphs, orders) =', per_size)
    print('Graph/shear symplectic and closed-neighbourhood checks:', graphs)
    print('Per-order reversed vertex-separation identities and valid bags:', orders)
    print('Independent introduce/forget pathwidth optima checked:', graphs)
    print('Source Corollary 1: exact rank/count cuts, all graph shears up to 4 modes:', rank_cuts)

    # Source Section V chain: sequential nearest-neighbour two-mode squeezers.
    # cosh(r)=5/3, sinh(r)=4/3 are exact rational values with cosh^2-sinh^2=1.
    chain_counts = []
    for n in range(2, 9):
        s = [[Fraction(int(i == j)) for j in range(2 * n)] for i in range(2 * n)]
        for v in range(n - 1):
            g = [[Fraction(int(i == j)) for j in range(2 * n)] for i in range(2 * n)]
            for r in range(2):
                i, j = 2 * v + r, 2 * (v + 1) + r
                g[i][i] = g[j][j] = Fraction(5, 3)
                g[i][j] = g[j][i] = Fraction(4 if r == 0 else -4, 3)
            s = matmul(g, s)
        assert matmul(matmul(transpose(s), omega(n)), s) == omega(n)
        expected_supports = [(1 << min(v + 2, n)) - 1 for v in range(n)]
        assert supports(s) == expected_supports
        forward, backward = tuple(range(n)), tuple(reversed(range(n)))
        f_profile = source_profile(supports(s), forward)
        b_profile = source_profile(supports(s), backward)
        assert f_profile == [2] * (n - 1) + [1]
        assert b_profile == list(range(n, 0, -1))
        fidelity_rank(s, forward)
        fidelity_rank(s, backward)
        chain_counts.append((n, max(f_profile), max(b_profile)))
        if n == 5:
            print('Source five-mode TMS-chain occupation profiles:', f_profile, b_profile)
    print('Source TMS-chain fidelity (modes, forward memory, reverse memory):', chain_counts)

    # Exact non-shear Gaussian matrices: local Fourier rotations and integral shears.
    additional = 0
    for n in range(1, 5):
        for mask in range(min(8, 1 << (n * (n - 1) // 2))):
            a = adjacency(n, mask)
            low = shear(a)
            fourier = omega(n)
            upper = transpose(shear([[2 * x for x in row] for row in a]))
            s = matmul(matmul(low, fourier), upper)
            assert matmul(matmul(transpose(s), omega(n)), s) == omega(n)
            for order in permutations(range(n)):
                additional += fidelity_rank(s, order)
    print('Source Corollary 1: additional exact non-shear rank/count cuts:', additional)
    assert graphs == 1099 and orders == 124469
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
