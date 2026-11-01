from math import isqrt

SOURCE = [0, 1, 2, 3, 3, 3, 5, 5, 4, 5, 5, 7, 7, 5, 7, 7, 7, 7, 8, 9, 7]


def crosses(n, x, y):
    corners = (x*x + y*y, (x+1)*(x+1) + y*y, x*x + (y+1)*(y+1), (x+1)*(x+1) + (y+1)*(y+1))
    return min(corners) < n < max(corners)

def quadrant_crosses(n, x, y):
    return x*x + y*y < n < (x+1)*(x+1) + (y+1)*(y+1)


def full_count(n):
    # A crossing square has a lower-left coordinate in this finite box:
    # if its lower-left norm is < n then |x|,|y| < sqrt(n).
    r = isqrt(n) + 2
    total = 0
    for x in range(-r, r + 1):
        for y in range(-r, r + 1):
            total += crosses(n, x, y)
    return total


def a(n):
    total = full_count(n)
    assert total % 4 == 0
    return total // 4


def quadrant_count(n):
    r = isqrt(n) + 2
    return sum(quadrant_crosses(n, x, y) for x in range(r + 1) for y in range(r + 1))


def is_sum_two_squares(n):
    r = isqrt(n)
    return any(x*x + y*y == n for x in range(r + 1) for y in range(r + 1))


# Fidelity against the displayed initial terms of OEIS A393562.
computed = [a(n) for n in range(len(SOURCE))]
assert computed == SOURCE, (computed, SOURCE)

# The four-axis reflection orbits reduce the count to the nonnegative quadrant.
for n in range(0, 80):
    assert full_count(n) == 4 * quadrant_count(n), n

# Expected settlement: every odd k > 1 has two distinct squared radii with
# the same value.  The no-sum-of-two-squares lemma gives Q(N) = 2s - 1
# whenever (s-1)^2 < N < s^2 and N is not a sum of two squares.
def witnesses(s):
    if s == 2:
        return 3, 4
    if s == 3:
        return 6, 7
    if s >= 4 and s % 2 == 0:
        return s*s - 1, s*s - 5
    if s >= 5 and s % 2 == 1:
        return s*s - 2, s*s - 6
    raise AssertionError(s)

for s in range(2, 251):
    u, v = witnesses(s)
    assert u != v and u > 0 and v > 0
    if s >= 3:
        assert (s - 1)**2 < u < s*s
        assert (s - 1)**2 < v < s*s
    if s >= 3:
        assert not is_sum_two_squares(u)
        assert not is_sum_two_squares(v)
    assert quadrant_count(u) == 2*s - 1
    assert quadrant_count(v) == 2*s - 1
    assert quadrant_count(u) == quadrant_count(v) == 2*s - 1

# Explicitly check the two small base pairs from the proposed argument.
assert (a(3), a(4)) == (3, 3)
assert (a(6), a(7)) == (5, 5)

print(f"fidelity_terms={len(SOURCE)}")
print("witnesses_checked_s=2..250")
print("ALL_CHECKS_OK")
