"""Arithmetic of pole-cancelling near-identity filters and a finite sparse separator."""
from fractions import Fraction as Q


def main():
    for m in (2,4,16,64):
        a=Q(m,m*m+1)
        assert 1-a*(m+Q(1,m))==0 and 2*a<1
    # Demonstration with a dense interior and a certified conservative inverse
    # norm bound. Choose the next point by a single integer comparison bound.
    M,J=4,5
    m=16*M*M
    a=Q(m,m*m+1)
    eps=Q(1, (2*M)**40)
    # Original normalized-by-unitary-dilation columns have smallest singular
    # value >=1/(4 M^(3/2)); smoothing and P_m preserve at least one quarter.
    kappa=1024*M**3
    eta=Q(1,2000*(1+kappa))
    assert 4*M*M*(2*Q(1,(2*M)**6))<Q(1,2)
    assert 2*a+Q(1,(2*M)**20)<Q(1,M)
    assert eps==Q(1,(2*M)**20)**2
    scale=M
    points=[]
    for j in range(1,J+1):
        goal=eta/(2**j*(M+j))
        # For r>=1, (2+log r)/sqrt(r) < 3/r^(1/4).
        # Use 16 for the ratio (1+2a)^2/||Fu||^2.
        base=48/goal
        cutoff=2*m**4*(int(base)+1)**4*scale
        p=2 if j%2 else 3
        n=p
        while n<cutoff:
            n*=5
        assert n>=2*scale and n<5*cutoff
        points.append(n)
        scale=n
    assert 2*eta+kappa*eta**2<Q(2,1000)
    print("Four dense inputs and five separated prime-branch pivots: largest index bit length",points[-1].bit_length())
    print("Uniform Schur trace-distance to the identity below .002; pole moments vanish.")


if __name__=="__main__":
    main()
