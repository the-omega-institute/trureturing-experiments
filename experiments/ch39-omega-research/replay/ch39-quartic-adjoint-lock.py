"""Small irreducibles, adjoint projector, and gauge masses in the quartic model."""
from fractions import Fraction as F
from itertools import combinations, product
from math import prod


def dimension(labels):
    n = len(labels)+1
    return prod(F(j-i+sum(labels[i:j]),j-i) for i in range(n) for j in range(i+1,n))


def small(n):
    # A nonnegative Dynkin label with too large a one-coordinate dimension can only grow.
    axes = []
    for i in range(n-1):
        vals = [0]*(n-1)
        while dimension(vals) < 24:
            vals[i] += 1
        axes.append(range(vals[i]))
    return [(tuple(z),int(dimension(z)),sum((j+1)*x for j,x in enumerate(z))%n)
            for z in product(*axes) if dimension(z)<24]


def verify():
    five,ten = small(5),small(10)
    assert sorted(d for _,d,_ in five)==[1,5,5,10,10,15,15]
    assert sorted(d for _,d,_ in ten)==[1,10,10]
    assert all(d==1 or k%5 for _,d,k in five+ten)
    assert dimension((1,0,0,1))==24
    c=(1,1,1,0,0)
    spec=[c[i]+c[j] for i,j in combinations(range(5),2)]
    assert [spec.count(t) for t in (2,1,0)] == [3,6,1]
    # Center of H is tested at (zeta_5 I_5, zeta_5^2 I_10).
    # A unit adjoint gap affects only the 12 color/weak generators.
    assert sum((c[i]-c[j])**2 for i,j in product(range(5),repeat=2)) == 12
    print("One real SU(5) adjoint, dimension 24, suffices and saturates the lower bound.")
    print("Exterior-square projector eigenvalue multiplicities 3,6,1.")


if __name__ == "__main__":
    verify()
