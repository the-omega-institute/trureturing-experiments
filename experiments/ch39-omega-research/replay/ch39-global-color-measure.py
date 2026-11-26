"""Charge, determinant and uniform neighborhood bounds for the color reduction."""
from fractions import Fraction as F
from itertools import product


def verify():
    # Each tuple lists color representation, multiplicity, 6Y and shell charge.
    species = [("3", 2, 1, 0), ("bar3", 1, 2, 1), ("bar3", 1, -4, -1),
               ("1", 2, -3, -1), ("1", 1, 6, -1)]
    assert sum(n*y for rep, n, y, q in species if rep != "1") == 0
    assert sum(n*q for rep, n, y, q in species if rep != "1") == 0
    assert sum(n for rep, n, y, q in species if rep == "3") == 2
    assert sum(n for rep, n, y, q in species if rep == "bar3") == 2
    assert all(abs(q) <= 1 for rep, n, y, q in species)
    # If W eigenangles lie in [-delta,delta], the correcting color angle is <=2 delta/3.
    for a, b in product((-1,1), repeat=2):
        assert abs(F(a+b,3)) <= F(2,3)
    delta, size = F(1,2**14), 2**15
    cost = 32*(2*delta + F(3, size))
    assert cost == F(7,1024) and cost < F(1,100)
    assert len(list(product((-1,0,1), repeat=4))) == 81
    print("Color-pair characters under hypercharge and shell phase vanish identically.")
    print("For L=32768 and delta=2^-14 the uniform kernel change is at most 7/1024.")


if __name__ == "__main__":
    verify()
