"""Relative spin eta phase on S^2 times an oriented lens space, in turns."""
from fractions import Fraction as F
from itertools import product


def square_half(a, b):
    # Inverse of the intersection matrix of the negative even A2 plumbing.
    return -F(a*a+a*b+b*b, 3)


def boundary_exponent(fluxes, discrete):
    # Three B copies carry E, three T copies carry E dual.
    reps = [(1, 1)]*3 + [(-1, 0), (-1, 1), (-1, -1)]
    total = sum(f*q*square_half(k*discrete, 0)
                for f in fluxes for q,k in reps)
    return total % 1


def verify():
    for a, b, i, j in product(range(-2, 3), repeat=4):
        assert (square_half(a-2*i+j, b+i-2*j)-square_half(a, b)).denominator == 1
    # A generator of the boundary cokernel has order three.
    assert (-2*(-2)+(-1), (-2)-2*(-1)) == (3, 0)
    for a in range(-5, 6):
        for b in range(-5, 6):
            assert (square_half(a, b)-square_half(a-b, 0)).denominator == 1
    for fluxes in ([1]+[0]*9, [1]*10, [4, -1]+[0]*8, [-2, 1]+[0]*8):
        for k in range(-4, 5):
            assert boundary_exponent(fluxes,k) == -F(sum(fluxes)*k*k, 3) % 1
    # Finite cyclic bar chains. The Bockstein is the carry in addition mod three.
    # Check the cocycle and its value on the normalized cycle sum_j [1|j|1].
    def cocycle(a, b, c):
        return a*((b+c)//3)
    for a,b,c,d in product(range(3), repeat=4):
        diff = (cocycle(b,c,d)-cocycle((a+b)%3,c,d)+cocycle(a,(b+c)%3,d)
                -cocycle(a,b,(c+d)%3)+cocycle(a,b,c))
        assert diff % 3 == 0
    boundary = {(a,b): 0 for a,b in product(range(1,3),repeat=2)}
    for j in range(3):
        for pair, sign in [((j,1),1),(((1+j)%3,1),-1),((1,(j+1)%3),1),((1,j),-1)]:
            if 0 not in pair:
                boundary[pair] += sign
    assert not any(boundary.values()) and sum(cocycle(1,j,1) for j in range(3)) == 1
    print("Spin filling half-square -1/3; faithful unit flux and covering unit flux",
          "both give 2/3 of a turn modulo one. Finite bar pairing is a generator.")


if __name__ == "__main__":
    verify()
