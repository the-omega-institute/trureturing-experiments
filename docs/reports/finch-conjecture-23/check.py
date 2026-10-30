#!/usr/bin/env python3
"""Exact checks for Finch, arXiv:2610.34300v1, Conjecture 23.

An order-k state collision represents the entire infinite tail. No assertion
about an infinite tail is inferred merely from a long scalar prefix.
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class Orbit:
    k: int
    values: tuple[int, ...]  # index zero is an unused sentinel
    first_endpoint: int
    repeat_endpoint: int

    @property
    def period(self) -> int:
        return self.repeat_endpoint - self.first_endpoint

    @property
    def tail_start(self) -> int:
        return self.first_endpoint - self.k + 1

    def value(self, n: int) -> int:
        assert n >= 1
        if n < self.tail_start:
            return self.values[n]
        return self.values[self.tail_start + (n - self.tail_start) % self.period]


def orbit(k: int) -> Orbit:
    assert k >= 2
    x = [0] + [0] * (k - 1) + [1]
    seen = {tuple(x[1:]): k}
    for n in range(k + 1, 200001):
        x.append((x[n - 1] + x[n - k] ** 2) % 3)
        state = tuple(x[n - k + 1:n + 1])
        if state in seen:
            j = seen[state]
            assert tuple(x[j - k + 1:j + 1]) == state
            result = Orbit(k, tuple(x), j, n)
            # Full recurrence fidelity over the stored certificate.
            assert all(x[t] == (x[t - 1] + x[t - k] ** 2) % 3
                       for t in range(k + 1, n + 1))
            # Every possible smaller period fails on a recurring cycle phase.
            assert all(any(result.value(result.tail_start + i + p)
                           != result.value(result.tail_start + i)
                           for i in range(result.period))
                       for p in range(1, result.period))
            return result
        seen[state] = n
    raise AssertionError(f"orbit cap exceeded at k={k}")


def preperiod(o: Orbit) -> int:
    # Beyond tail_start all phases recur forever; before it there are finitely
    # many possible disagreements. The strict source inequality is n > mu.
    mismatches = [n for n in range(1, o.tail_start)
                  if o.value(n + o.period) != o.value(n)]
    mu = max(mismatches, default=0)
    assert all(o.value(n + o.period) == o.value(n)
               for n in range(mu + 1, o.repeat_endpoint + o.period + 1))
    if mu:
        assert o.value(mu + o.period) != o.value(mu)
    return mu


def capture(o: Orbit) -> int:
    assert o.period == 4
    gaps = (o.k - 1 + 3) // 4
    assert gaps >= 1
    # An eventual four-cycle containing a zero supplies such a run within a
    # complete period after max(tail_start, k+1), if there was no earlier run.
    bound = max(o.tail_start, o.k + 1) + o.period + 4 * gaps
    for z in range(o.k + 1 + 4 * gaps, bound + 1):
        r = z - 4 * gaps
        if all((o.value(r + i) == 0) == (i % 4 == 0)
               for i in range(4 * gaps + 1)):
            return z
    raise AssertionError(f"no source capture run at k={o.k}")


def main() -> None:
    # Fidelity: reproduce every period in the main-text recurrence (A) table.
    expected = {2: 4, 3: 4, 4: 9, 5: 19, 6: 4, 7: 4,
                8: 22, 9: 36, 10: 4, 11: 4, 12: 45}
    for k, period in expected.items():
        actual = orbit(k).period
        assert actual == period, (k, actual, period)
    print("FIDELITY main-text N_A(k), k=2..12:",
          ",".join(str(expected[k]) for k in range(2, 13)))
    three = orbit(3)
    assert preperiod(three) == 5
    assert capture(three) == 11
    print("FIDELITY Addendum evidence: k=3, mu=5, z=11, z-mu=6")

    # Independently test the proposed closed four-cycle in every transition
    # phase, then compare against the recurrence's separately generated orbit.
    block = (0, 1, 1, 2)
    value = lambda n: block[(n - 1) % 4]
    assert value(1) == 0 and value(2) == 1
    assert all(value(n) == (value(n - 1) + value(n - 2) ** 2) % 3
               for n in range(3, 7))
    two = orbit(2)
    assert (two.first_endpoint, two.repeat_endpoint) == (2, 6)
    assert all(two.value(n) == value(n) for n in range(1, 81))
    assert all(any(value(n + p) != value(n) for n in range(1, 5))
               for p in (1, 2, 3))
    prefix = tuple(two.value(n) for n in range(1, 10))
    assert prefix == (0, 1, 1, 2, 0, 1, 1, 2, 0)
    assert [n for n in range(3, 10) if two.value(n) == 0] == [5, 9]
    assert preperiod(two) == 0 and capture(two) == 9
    assert capture(two) - preperiod(two) == 9 != 6
    print("COUNTEREXAMPLE k=2: x_1..x_9=" + ",".join(map(str, prefix)))
    print("CERTIFICATE state (x_1,x_2)=(x_5,x_6)=(0,1); least period=4")
    print("COUNTEREXAMPLE k=2: G=1, first admissible zeros=5,9; mu=0; z=9; z-mu=9 != 6")

    # Additional source-listed evidence, without asserting a corrected uniform
    # conjecture or an infinite exceptional family.
    for k in (3, 6, 7, 10, 11, 14, 15, 19):
        o = orbit(k)
        assert o.period == 4
        mu, z = preperiod(o), capture(o)
        assert z - mu == 6
        print(f"SOURCE_LISTED k={k}: N_A=4, mu={mu}, z={z}, z-mu={z-mu}")
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
