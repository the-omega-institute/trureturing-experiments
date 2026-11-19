"""Exact checks for the Cambie–Jooken–Wagner Question 1.9 candidate."""

from fractions import Fraction
from itertools import combinations, product
from math import comb


# Polynomials in ascending order of powers of t.
NI = (12, -19, 256, 403, 166)
RI3 = (60, 115, 542, 6531, 10027, 4190)
NT = (9, 31, 471, 627, 232)
RT3 = (54, 242, 3748, 14388, 16505, 5854)
D = (-108, 1017, -17499, -17792, -356655, -1199435,
     -1461298, -762936, -145598, 316)


def evaluate(coefficients, t):
    value = 0
    for coefficient in reversed(coefficients):
        value = value * t + coefficient
    return value


def multiply(left, right):
    coefficients = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            coefficients[i + j] += a * b
    return coefficients


def literal_counts(size, edges):
    """Enumerate vertex subsets for I and actual edge subsets for T."""
    neighbors = [0] * size
    endpoints = []
    for u, v in edges:
        neighbors[u] |= 1 << v
        neighbors[v] |= 1 << u
        endpoints.append((1 << u) | (1 << v))

    count_i = order_i = 0
    for vertices in range(1, 1 << size):
        pending = vertices & -vertices
        reached = pending
        while pending:
            vertex = pending & -pending
            fresh = neighbors[vertex.bit_length() - 1] & vertices & ~reached
            reached |= fresh
            pending ^= vertex
            pending |= fresh
        if reached == vertices:
            count_i += 1
            order_i += vertices.bit_count()

    count_t = order_t = size  # The edgeless trees are exactly the singletons.
    for edge_mask in range(1, 1 << len(edges)):
        chosen = [edge for i, edge in enumerate(endpoints) if edge_mask & (1 << i)]
        vertices = 0
        for edge in chosen:
            vertices |= edge
        order = vertices.bit_count()
        if len(chosen) != order - 1:
            continue
        reached = chosen.pop()
        while chosen:
            rest = []
            for edge in chosen:
                if edge & reached:
                    reached |= edge
                else:
                    rest.append(edge)
            if len(rest) == len(chosen):
                break
            chosen = rest
        if reached == vertices:
            count_t += 1
            order_t += order
    return count_i, order_i, count_t, order_t


def graph(t):
    m = 10 * t
    edges = [(i, (i + 1) % m) for i in range(m)]
    for anchor, start in ((0, m), (2 * t, m + t)):
        edges.append((anchor, start))
        edges.extend((i, i + 1) for i in range(start, start + t - 1))
    x = m + 2 * t
    edges += [(t, x), (x, x + 1), (x + 1, x + 2), (x + 2, x)]
    return x + 3, edges


def cycle_moments(m, marks, trees):
    if not marks:
        if trees:
            return m * m, m * m * (m + 1) // 2
        return m * (m - 1) + 1, m * m * (m - 1) // 2 + m
    gaps = [v - u for u, v in zip(marks, marks[1:])]
    gaps.append(m - marks[-1] + marks[0])
    n = (m if trees else 1) + sum(comb(d, 2) for d in gaps)
    return n, m * n - sum(comb(d + 1, 3) for d in gaps)


def marked_gap_counts(t, trees):
    """Evaluate the gap formula and the attachment value/derivative product."""
    triangle_count, triangle_order = (9, 18) if trees else (7, 12)
    count = t * (t + 1) + triangle_count
    order = t * (t + 1) * (t + 2) // 3 + triangle_order
    # Values and first derivatives of p_t - 1 and h - 1 at w = 1.
    path = (t, t * (t + 1) // 2)
    triangle = (6, 14) if trees else (4, 8)
    attachments = ((0, path), (t, triangle), (2 * t, path))
    for subset_size in range(4):
        for subset in combinations(attachments, subset_size):
            marks = [vertex for vertex, _ in subset]
            cycle_count, cycle_order = cycle_moments(10 * t, marks, trees)
            weight, derivative = 1, 0
            for _, (value, slope) in subset:
                derivative = derivative * value + weight * slope
                weight *= value
            count += cycle_count * weight
            order += cycle_order * weight + cycle_count * derivative
    return count, order


def direct_interval_counts(t, trees):
    """Count cycle intervals by their actual marks, without gap formulas."""
    m = 10 * t
    triangle = [0, 1, 2, 2, 3, 3, 3] if trees else [0, 1, 2, 2, 3]
    no_cycle = 2 * [k for k in range(1, t + 1) for _ in range(t - k + 1)]
    no_cycle += [1, 1, 1, 2, 2, 2] + ([3, 3, 3] if trees else [3])
    count, order = len(no_cycle), sum(no_cycle)

    intervals = [(start, length, 1) for start in range(m) for length in range(1, m)]
    intervals.append((0, m, m if trees else 1))
    for start, length, multiplicity in intervals:
        prefixes = [range(t + 1) if (anchor - start) % m < length else (0,)
                    for anchor in (0, 2 * t)]
        rooted_orders = triangle if (t - start) % m < length else [0]
        for extra in product(*prefixes, rooted_orders):
            count += multiplicity
            order += multiplicity * (length + sum(extra))
    return count, order


def main():
    for n in range(1, 9):
        ni, ri, nt, rt = literal_counts(n, [(i, i + 1) for i in range(n - 1)])
        assert Fraction(rt, nt) == Fraction(ri, ni) == Fraction(n + 2, 3)
    print("Source fidelity: Theorem 1.2 path equality and the tree-case equality on P_1,...,P_8")

    literal = literal_counts(*graph(1))
    assert literal == (*direct_interval_counts(1, False), *direct_interval_counts(1, True))
    print(f"Literal G_1 (N_I, R_I, N_T, R_T): {literal}")

    for t in range(1, 5):
        for trees in (False, True):
            assert marked_gap_counts(t, trees) == direct_interval_counts(t, trees)
    print("Direct cycle-interval enumeration agrees with the marked-gap formula for t=1,2,3,4")

    # The gap formulas have degrees at most 4 for N and 5 for 3R:
    # cycle N and R have degrees at most 2 and 3, and attachments
    # multiply by at most two path factors, each of degree 1 in N.
    # A derivative can replace one path factor by a degree-2 factor.
    for t in range(1, 7):
        for trees, n_poly, r3_poly in ((False, NI, RI3), (True, NT, RT3)):
            count, order = marked_gap_counts(t, trees)
            assert count == evaluate(n_poly, t)
            assert 3 * order == evaluate(r3_poly, t)
    assert D == tuple(a - b for a, b in zip(multiply(RI3, NT), multiply(RT3, NI)))
    print("Four polynomial identities checked at six distinct t; cross-product verified coefficientwise")

    shifted = tuple(sum(D[j] * comb(j, i) * 500 ** (j - i) for j in range(i, len(D)))
                    for i in range(len(D)))
    assert all(coefficient > 0 for coefficient in shifted)
    print(f"D(500+u) coefficients, constant first: {shifted}")

    ni, ri = marked_gap_counts(500, False)
    nt, rt = marked_gap_counts(500, True)
    difference = Fraction(ri, ni) - Fraction(rt, nt)
    assert difference == Fraction(1769250196421378505239933, 18998398362554330511481326)
    assert difference > 0
    assert Fraction(D[-1], 3 * NI[-1] * NT[-1]) == Fraction(79, 28884)
    print(f"G_500: av - mu = {difference}")
    print("Linear coefficient of av - mu: 79/28884")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
