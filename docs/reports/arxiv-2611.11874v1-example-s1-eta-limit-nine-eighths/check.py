"""Exact checks of Corenflos, arXiv:2611.11874v1, Example S1."""

from fractions import Fraction as F
from itertools import product
from math import comb, lcm

NULL = (3, 4, 5, 6, 7)
CAP = F(9, 8)
LAST = 200


def mass(a, t, k):
    """Integer numerator of the probability of a specified word over 10**t."""
    return a**k * (10 - a)**(t - k)


def reward(t, k):
    return F(mass(2, t, k) + mass(8, t, k), mass(3, t, k) + mass(7, t, k))


def explicit_policy_values(prefix, horizon, theta):
    """All possible conditional values on the full ordered prefix tree."""
    t = len(prefix)
    stop = reward(t, sum(prefix))
    if t == horizon:
        return {stop}
    left = explicit_policy_values(prefix + (0,), horizon, theta)
    right = explicit_policy_values(prefix + (1,), horizon, theta)
    return {stop} | {(1 - theta) * a + theta * b for a, b in product(left, right)}


def exact_snell(horizon, stops, common):
    """Backward induction with a fixed integer scale at each tree depth."""
    values = {a: stops[horizon] for a in NULL}
    for t in range(horizon - 1, -1, -1):
        values = {
            a: [max(stop, (10 - a) * v[k] + a * v[k + 1])
                for k, stop in enumerate(stops[t])]
            for a, v in values.items()
        }
    return {a: F(v[0], common * 10**LAST) for a, v in values.items()}


def first_hit_check(rows, integer_rows):
    """Unnormalized P(0.7) prefix masses of paths still above -1."""
    survivors = [1]
    stopped = F(0)
    for t, row in enumerate(rows):
        survivors = (survivors + [0] * (len(row) - len(survivors)))
        if t % 2:
            k = t // 2
            stopped += F(survivors[k], 10**t) * row[k]
            survivors[k] = 0
        den, values = integer_rows[t]
        yield stopped + F(sum(w * f for w, f in zip(survivors, values)), 10**t * den)
        next_mass = [0] * (t + 2)
        for k, w in enumerate(survivors):
            next_mass[k] += 3 * w
            next_mass[k + 1] += 7 * w
        survivors = next_mass


def main():
    rows = [[reward(t, k) for k in range(t + 1)] for t in range(LAST + 1)]
    integer_rows = []
    for t, row in enumerate(rows):
        den = lcm(*(f.denominator for f in row))
        values = [f.numerator * (den // f.denominator) for f in row]
        weights = [comb(t, k) * f for k, f in enumerate(values)]
        unit = den * 10**t
        for a in NULL:
            mean_numerator = sum(mass(a, t, k) * f for k, f in enumerate(weights))
            assert mean_numerator <= unit, ("numeraire feasibility", t, a)
            if a in (3, 7) or t <= 1:
                assert mean_numerator == unit, (t, a)
            else:
                assert mean_numerator < unit, (t, a)
        integer_rows.append((den, values))
    common = lcm(*(den for den, _ in integer_rows))
    # One depth earlier has 10 times the denominator, absorbing Bernoulli /10.
    stops = []
    for t, (den, values) in enumerate(integer_rows):
        scale = (common // den) * 10**(LAST - t)
        stops.append([v * scale for v in values])
    print("Numeraire: all five null constraints verified exactly for t=0..200.")

    fidelity = F(3, 10) * rows[1][0] + F(7, 10) * (
        F(3, 10) * rows[2][1] + F(7, 10) * rows[2][2]
    )
    assert fidelity == F(30, 29)
    print("Example S1 stopping rule at P(0.7): E[E*_{1+1(X1=1)}] = 30/29.")

    previous = F(0)
    for t, hit in enumerate(first_hit_check(rows, integer_rows)):
        values = exact_snell(t, stops, common)
        if t <= 4:
            for a in NULL:
                assert values[a] == max(explicit_policy_values((), t, F(a, 10)))
        eta = max(values.values())
        assert values[3] == values[7] == eta, ("endpoint maximizers", t)
        assert previous <= eta, ("monotonicity", t)
        lower = CAP - F(17, 21)**t - F(1, 8) * F(4, 5)**t
        assert lower <= hit <= values[7] <= CAP, ("envelope", t)
        if t <= 4:
            print(f"eta_{t} = {eta}")
        if t in (50, 200):
            lo, hi = {50: (F(465, 10**9), F(475, 10**9)),
                      200: (F(325, 10**24), F(335, 10**24))}[t]
            gap = CAP - eta
            assert lo < gap < hi, ("source rounding", t)
            print(f"Source gap at T={t}: {lo} < 9/8 - eta_{t} < {hi}.")
        previous = eta
    print("Full ordered-tree policy enumeration agrees for T=0..4 at all five nulls.")
    print("Exact Snell values, endpoint attainment, monotonicity and both bounds verified for T=0..200.")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
