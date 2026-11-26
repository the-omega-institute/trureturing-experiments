"""Uniform-in-kinetic-strength Haar caps and the finite low-mode count."""
from fractions import Fraction
from math import comb


def certify():
    dimensions = (9, 7, 5, 3)
    assert sum(dimensions) == 24 and sum(dimensions[1:]) == 15
    for dimension in dimensions:
        assert 2**dimension >= 2*(dimension+1)
    assert 4**2*3 < 8**2
    assert Fraction(19, 10)*(8*Fraction(7, 4)/56 + Fraction(16, 64)) < 1
    assert (Fraction(19, 10)-1)**2 > Fraction(4, 5) and Fraction(7, 4)**2 > 3
    # Adaptive, correlated block successes with conditional probability at least p.
    p = Fraction(1, 7)
    law = {0: Fraction(1)}
    for size in range(1, 9):
        updated = {}
        for count, mass in law.items():
            prob = p + (1-p)*Fraction((size+count) % 4, 4)
            updated[count] = updated.get(count, 0) + (1-prob)*mass
            updated[count+1] = updated.get(count+1, 0) + prob*mass
        law = updated
        assert sum(law.values()) == 1
        for threshold in range(1, size+1):
            lower_tail = sum(comb(size,j)*p**j*(1-p)**(size-j) for j in range(threshold))
            assert sum(mass for j, mass in law.items() if j < threshold) <= lower_tail
    eta = Fraction(1, 100)
    radius, side = 5600, 11203
    assert 56/eta == radius and 2*radius+3 == side
    assert 32/(eta/64) < 2**18 and 23*14+24*18 == 754
    assert 4*754 == 3016 and 2**4 >= 10
    print("Verified 24 cap powers, both spin-color counts' common tolerance and box budget.")


if __name__ == "__main__":
    certify()
