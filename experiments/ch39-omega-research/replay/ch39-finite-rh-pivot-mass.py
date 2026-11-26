"""Six dyadic RH Gram readings determine three threshold comparisons."""
from fractions import Fraction as Q
from itertools import permutations

bins = [
    [],
    [79164649776018547114746344827],
    [78976713806822045219190780366,79165616124469926724257714927],
    [78671837523006435877920812156,78982331081163686403340645257,79167467220854044217387473963],
]
den = 2**96


def read(i,j):
    if i==j:
        return Q(1),Q(1)
    k = bins[max(i,j)][min(i,j)]
    return Q(k,den),Q(k+1,den)


def determinant(indices):
    lo = hi = Q(0)
    n = len(indices)
    for p in permutations(range(n)):
        a = b = Q(1)
        for i,j in enumerate(p):
            lower,upper = read(indices[i],indices[j])
            a *= lower
            b *= upper
        if sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))%2:
            lo -= b
            hi -= a
        else:
            lo += a
            hi += b
    return lo,hi


def main():
    target = [(Q("24.97936"),Q("24.97938")),(Q("265.97205"),Q("265.97207")),
              (Q("2609.30190"),Q("2609.30192"))]
    prev = determinant(range(1))
    for n,(lo,hi) in enumerate(target,2):
        curr = determinant(range(n))
        assert curr[0]>0 and lo**2<prev[0]/curr[1] and prev[1]/curr[0]<hi**2
        prev = curr
    assert target[1][1]<Q("266.6444299")<target[2][0]
    assert target[2][0]**2*8*14006757711**2>50017863218932**2
    all_three = determinant(range(3))
    for omit,lo,hi in [(0,"263.94171","263.94173"),(1,"528.89869680","528.89869682")]:
        minor = determinant([i for i in range(3) if i!=omit])
        assert Q(lo)**2<minor[0]/all_three[1] and minor[1]/all_three[0]<Q(hi)**2
    print("Inverse square-root Schur pivots:",[(str(x),str(y)) for x,y in target])
    print("The third node misses the central threshold; the fourth decouples both charged thresholds.")
    print("The same three nodes give inverse residual norm between 528.89869680 and 528.89869682 at the middle node.")


if __name__=="__main__":
    main()
