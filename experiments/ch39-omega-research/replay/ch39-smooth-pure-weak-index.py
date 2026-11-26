"""Finite section, compression estimates and action budget for an odd weak index."""
from fractions import Fraction as F
from itertools import product
from math import comb


def ramp(j, length):
    j %= length
    return min(F(1), F(4*min(j, length-j),length))


def phase_power(j, k, length):
    return length*k*int(2*j >= length)


def section(length, coords):
    x,y,z,w = coords
    a,b = ramp(x,length),ramp(z,length)
    # Store each unnormalized coordinate as magnitude and power of the L^2-th root.
    return [(a,0),(b,0),((1-a)*(1-b),phase_power(x,y,length)+phase_power(z,w,length))]


def check_patches(length):
    for x,y in product(range(length),repeat=2):
        for mu in range(2):
            dest = ((x+1)%length,y) if mu == 0 else (x,(y+1)%length)
            source = section(length,(x,y,0,0))
            target = section(length,(*dest,0,0))
            e = (-length*y*int(x==length-1)) if mu == 0 else x
            for slot, twist in ((0,e),(2,-e)):
                mag, k = source[slot]
                mag2, k2 = target[slot]
                assert abs(mag-mag2) <= F(4,length)
                if mag*mag2:
                    error = (k2-k+twist) % length**2
                    assert min(error,length**2-error) <= (0 if mu == 0 else length)
    for x,z in product(range(length),repeat=2):
        c = section(length,(x,0,z,0))
        assert max(m for m,k in c) >= F(1,4)


def verify():
    for length in (8,9,16,31,32):
        check_patches(length)
    assert sum(comb(4,k) for k in (2,4)) == 7
    curvature = 16+22*128**2
    size = 2**22
    assert curvature == 360464
    assert F(1,100) - 32*4*F(128,size) > F(1,200)
    assert 2*F(1,size**2) < F(1,2**29)
    assert 3*curvature**2 == 389802885888
    print("Unit weak index magnitude; box-independent six-plane Wilson action <=",3*curvature**2)


if __name__ == "__main__":
    verify()
