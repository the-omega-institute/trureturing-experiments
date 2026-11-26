"""Exact partially transposed Gram spectra for charged weak singlet currents."""
from collections import Counter
from math import comb, factorial
from flint import fmpq, fmpq_mat, fmpq_poly


def spin_four():
    # Projection onto the normalized coefficient vectors of a pair of binary octics.
    n = 8
    assert all(sum(comb(n,i)*comb(n,m-i) for i in range(n+1) if 0 <= m-i <= n)
               == comb(2*n,m) for m in range(2*n+1))
    characteristic = fmpq_poly([1])
    for difference in range(-n,n+1):
        block = [(i,i-difference) for i in range(n+1) if 0 <= i-difference <= n]
        mat = fmpq_mat([[fmpq(comb(n,i)*comb(n,l),comb(2*n,i+j))
                         for k,j in block] for i,l in block])
        characteristic *= mat.charpoly()
    expected = fmpq_poly([1])
    for l in range(n+1):
        eigen = fmpq((2*n+1)*factorial(n)**2, factorial(n-l)*factorial(n+l+1))
        expected *= fmpq_poly([-eigen,1])**(2*l+1)
    assert expected == characteristic and eigen == fmpq(1,12870)
    print("17 orthonormal charged forms, partial transpose floor", eigen)


def tensor_frames():
    # For normalized symmetric-matrix frames in C^3 tensor C^3, the partial transpose
    # of P_sym is (I+|Omega><Omega|)/2. Omitting Omega/sqrt(3) subtracts swap/3.
    for left, right, rank, lower in [(5,5,25,fmpq(1,36)), (5,6,30,fmpq(1,12))]:
        spectra = {}
        for dimension in (5,6):
            part = Counter()
            for mult, swap, trace in ((5,1,0), (3,-1,0), (1,1,3)):
                value = fmpq(1+trace,2)-(fmpq(swap,3) if dimension==5 else 0)
                part[value] += mult
            spectra[dimension] = part
        assert left*right == rank
        assert min(spectra[left])*min(spectra[right]) == lower
        print(rank, "orthonormal charged forms, partial transpose floor", lower)


if __name__ == "__main__":
    spin_four()
    tensor_frames()
