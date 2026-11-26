"""Four finite backgrounds certify the sharp determinant-character spectator count."""
from collections import Counter
from itertools import product


def three_flux_rows(n, q):
    return tuple(n*int((n*k+q) % 3 != 0) for k in range(3))


def add_rows(fields):
    return tuple(sum(three_flux_rows(n,q)[k] for n,q in fields) for k in range(3))


def verify():
    original = [(1,-1), (1,-1), (-1,-1), (-1,0)]
    mirrors = [(-n, -q) for n,q in original]
    original_rows = add_rows(original)
    assert original_rows == (1,-2,1)
    forced = [(1,0), (1,1), (-2,-1)]
    assert add_rows(forced) == tuple(-x for x in original_rows)
    index5 = (0,50,25,25,50)
    assert sum(n*index5[n%5] for n,q in original) == 0
    assert sum(n*index5[n%5] for n,q in forced) == 50
    counts = Counter((n,q%3) for n,q in original+mirrors)
    assert all(count == counts[-n,(-q)%3] for (n,q),count in counts.items())
    for n,q in product(range(-8,9),range(3)):
        row = three_flux_rows(n,q)
        if n % 3:
            absent = (-q*pow(n,-1,3)) % 3
            assert row[absent] == 0 and sum(row) == 2*n
        else:
            assert len(set(row)) == 1 and row[0]%3 == 0
    print("Normalized central charges", original_rows, "and the forced triple fails the fifth-root flux.")


if __name__ == "__main__":
    verify()
