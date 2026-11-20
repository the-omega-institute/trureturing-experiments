"""Exact checks for the exponential sphere-functional at k = 1."""

from collections import Counter, deque
from fractions import Fraction
from math import factorial, isqrt

UNIT = 1 << 48
ONE = (UNIT, UNIT)
ZERO = (0, 0)


def ceil_div(n, d):
    return -(-n // d)


def enclose(lo, hi):
    return (
        UNIT * lo.numerator // lo.denominator,
        ceil_div(UNIT * hi.numerator, hi.denominator),
    )


def add(x, y):
    return x[0] + y[0], x[1] + y[1]


def subtract(x, y):
    return x[0] - y[1], x[1] - y[0]


def scale(x, n):
    if n < 0:
        return x[1] * n, x[0] * n
    return x[0] * n, x[1] * n


def multiply(x, y):
    assert x[0] >= 0 and y[0] >= 0
    return x[0] * y[0] // UNIT, ceil_div(x[1] * y[1], UNIT)


def divide(x, y):
    assert x[0] >= 0 and y[0] > 0
    return x[0] * UNIT // y[1], ceil_div(x[1] * UNIT, y[0])


def exp_minus_one():
    terms = [Fraction((-1) ** k, factorial(k)) for k in range(28)]
    # Alternating series: odd partial sum below, even partial sum above.
    return enclose(sum(terms), sum(terms[:-1]))


def log_near_one(value):
    # For 1 <= x <= 3, z = (x-1)/(x+1) <= 1/2 and
    # ln(x) = 2 * sum(z^(2k+1)/(2k+1), k >= 0).
    assert UNIT <= value[0] <= value[1] <= 3 * UNIT
    z = divide(subtract(value, ONE), add(value, ONE))
    z_squared = multiply(z, z)
    power = z
    lower = upper = 0
    for k in range(26):
        denominator = 2 * k + 1
        lower += power[0] // denominator
        upper += ceil_div(power[1], denominator)
        power = multiply(power, z_squared)
    # 2*(1/2)^53 / (53*(1-(1/2)^2)) < 1/UNIT.
    assert 8 * UNIT < 3 * 53 * (1 << 53)
    return 2 * lower, 2 * upper + 1


LOG_TWO = log_near_one(scale(ONE, 2))


def log_interval(value):
    assert value[0] > 0
    exponent = value[0].bit_length() - UNIT.bit_length()
    if exponent >= 0:
        divisor = 1 << exponent
        reduced = value[0] // divisor, ceil_div(value[1], divisor)
    else:
        reduced = scale(value, 1 << -exponent)
    return add(log_near_one(reduced), scale(LOG_TWO, exponent))


def power_tables(n, q):
    powers = [ONE]
    geometric = [ZERO]
    for _ in range(n):
        powers.append(multiply(powers[-1], q))
        geometric.append(add(geometric[-1], powers[-1]))
    return powers, geometric


def star_groups(length_counts, powers, geometric):
    center = ZERO
    for length, copies in length_counts:
        center = add(center, scale(geometric[length], copies))
    yield 1, center
    for length, copies in length_counts:
        # Exclude the current branch when summing through the center.
        others = subtract(center, geometric[length])
        others = max(others[0], 0), others[1]
        for depth in range(1, length + 1):
            same_arm = add(geometric[depth], geometric[length - depth])
            yield copies, add(same_arm, multiply(powers[depth], others))


def balanced(n, b):
    length, extra = divmod(n - 1, b)
    return [(length, b - extra)] + ([(length + 1, extra)] if extra else [])


def entropy(groups):
    groups = list(groups)
    total = ZERO
    for count, score in groups:
        total = add(total, scale(score, count))
    result = ZERO
    for count, score in groups:
        p = divide(score, total)
        negative_log = log_interval(divide(ONE, p))
        result = add(result, scale(multiply(p, negative_log), count))
    return result


def bfs_rows(graph):
    for start in range(len(graph)):
        distance = {start: 0}
        queue = deque([start])
        while queue:
            vertex = queue.popleft()
            for neighbor in graph[vertex]:
                if neighbor not in distance:
                    distance[neighbor] = distance[vertex] + 1
                    queue.append(neighbor)
        assert len(distance) == len(graph)
        del distance[start]
        yield Counter(distance.values())


def spider(length_counts):
    graph = [[]]
    for length, copies in length_counts:
        for _ in range(copies):
            previous = 0
            for _ in range(length):
                vertex = len(graph)
                graph[previous].append(vertex)
                graph.append([previous])
                previous = vertex
    return graph


def check_small_distance_sums():
    q = Fraction(1, 3)
    comparisons = 0
    for n in range(2, 17):
        exact_powers = [q ** j for j in range(n)]
        exact_geometric = [sum(exact_powers[1:j+1]) for j in range(n)]
        # Check at an exact rational q against an independent graph traversal.
        for length_counts in [balanced(n, b) for b in range(1, n)] + [
            [(1, m), (n-m-1, 1)] for m in range(1, n-2)
        ]:
            scores = sorted(
                sum(count * exact_powers[d] for d, count in row.items())
                for row in bfs_rows(spider(length_counts))
            )
            center = sum(count * exact_geometric[l] for l, count in length_counts)
            formula = [center]
            for length, copies in length_counts:
                for depth in range(1, length+1):
                    score = exact_geometric[depth] + exact_geometric[length-depth]
                    score += exact_powers[depth] * (center - exact_geometric[length])
                    formula.extend([score] * copies)
            assert scores == sorted(formula)
            comparisons += 1
    print(f"Exact breadth-first distance check: {comparisons} balanced stars and brooms")


def check_fidelity(q):
    n = 7
    cycle = [[(i-1) % n, (i+1) % n] for i in range(n)]
    rows = list(bfs_rows(cycle))
    assert all(row == {1: 2, 2: 2, 3: 2} for row in rows)
    # Equation (6c): every row has the identical nonzero score 6*(1+q+q^2).
    # Cancel its common symbolic factor to check the coefficient of log(7).
    probability = Fraction(1, n)
    assert sum(probability for _ in rows) == 1
    print("Source Proposition 4.1 on C_7: normalized weights 1/7, entropy log(7)")

    powers, geometric = power_tables(20, q)
    figure = balanced(20, 7)
    assert figure == [(2, 2), (3, 5)]
    figure_entropy = entropy(star_groups(figure, powers, geometric))
    for b in range(1, 20):
        if b != 7:
            assert figure_entropy[1] < entropy(star_groups(balanced(20, b), powers, geometric))[0]
    print("Source Figure 5: its order-20 tree uniquely minimizes among balanced spiders")


def check_family(n, q):
    powers, geometric = power_tables(n, q)
    log_n = log_interval(scale(ONE, n))
    # B*n*ln(n)/q^2 = 2*n*ln(n)/(q*(1-q)).
    square = divide(scale(log_n, 2*n), multiply(q, subtract(ONE, q)))
    m = isqrt(square[0] // UNIT)
    assert m == isqrt(square[1] // UNIT)
    tail = n - m - 1
    assert 2 <= m and 3 <= tail
    broom = entropy(star_groups([(1, m), (tail, 1)], powers, geometric))
    best_lower = min(
        entropy(star_groups(balanced(n, b), powers, geometric))[0]
        for b in range(1, n)
    )
    gap = Fraction(best_lower - broom[1], UNIT)
    assert gap > 0
    print(f"n={n}, m={m}, tail={tail}: all {n-1} arm counts, gap >= {gap} nats")


def main():
    q = exp_minus_one()
    check_small_distance_sums()
    check_fidelity(q)
    for n in (1024, 2048, 4096):
        check_family(n, q)
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
