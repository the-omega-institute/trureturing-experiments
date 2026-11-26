"""Finite character sums for weak-Cartan charges of a left-handed 5-bar plus 10."""
from itertools import product


def weak(index5, index10, a, b):
    n = len(index5)
    return index10[(a+b) % n]+index10[(a-b) % n]-2*index10[a]-index10[b]-index5[b]


def direct(index5, index10, charges):
    n = len(index5)
    weights = [-index5[-x % n]+sum(index10[(x+y) % n] for j, y in enumerate(charges) if i != j)
               for i, x in enumerate(charges)]
    return weights[3]-weights[4]


def even_tables(n):
    for values in product((-1, 0, 1), repeat=n//2):
        half = (0, *values)
        yield tuple(half[min(k, n-k)] for k in range(n))


def average_check():
    for n in range(2, 8):
        for f, g in product(even_tables(n), repeat=2):
            for b in range(n):
                row = [weak(f, g, a, b) for a in range(n)]
                assert all(row[a] == direct(f, g, (a, -a, b, -b, 0)) for a in range(n))
                assert row[0] == g[b]-f[b] and sum(row) == -n*(g[b]+f[b])
                assert max(map(abs, row)) >= abs(g[b])+abs(f[b])


def photon_check():
    indices = [0, 50, 25, 25, 50]
    for a in range(5):
        assert direct(indices, indices, (a, a, a, -3*a, 0)) == 3*indices[2*a % 5]-3*indices[a]-indices[3*a % 5]
    assert weak(indices, indices, 1, 1) == -175
    for n in range(2, 11):
        for f in even_tables(n):
            weak_zero = all(3*(f[2*a % n]-f[a]) == f[3*a % n] for a in range(n))
            exceptions = all(f[2*a % n] == f[a] and f[3*a % n] == 0 for a in range(n))
            assert weak_zero == exceptions


def bounded_charge_law():
    # Unwrapped charges permit a nonzero quadratic index and cancel each weak charge.
    indices = {i: 3*i*i for i in range(-20, 21)}
    for a, b in product(range(11), repeat=2):
        assert indices[a+b]+indices[a-b] == 2*indices[a]+2*indices[b]
    for charges in product(range(-2, 3), repeat=4):
        q = (*charges, -sum(charges))
        weights = [-indices[-x]+sum(indices[x+y] for j, y in enumerate(q) if j != i)
                   for i, x in enumerate(q)]
        assert set(weights) == {3*sum(x*x for x in q)}


if __name__ == "__main__":
    average_check()
    photon_check()
    bounded_charge_law()
    print("Verified the finite averages, electromagnetic exceptions and unwrapped quadratic law.")
