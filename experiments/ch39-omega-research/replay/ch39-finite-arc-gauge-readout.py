"""Spectral projectors of the six patch readings and an explicit even-occupation model."""
from collections import Counter
from itertools import combinations

from sympy import Matrix, Rational as Q, eye, ones


def main():
    labels = list(combinations(range(4),2))
    J = Matrix([[int(set(p).isdisjoint(q)) for q in labels] for p in labels])
    P1 = ones(6)/6
    P3 = (eye(6)-J)/2
    P2 = eye(6)-P1-P3
    M = Matrix([[Q(int(len(set(p)&set(q))==1),4) for q in labels] for p in labels])
    assert M == P1-P2/2
    assert all(P**2==P and P.rank()==k for k,P in ((1,P1),(3,P3),(2,P2)))
    assert M*P3 == 0*P3
    # Central circle with primitive integral weights -2 on 3 and 3 on 2.
    states = Counter()
    traces = [0,0,0,0]
    for k in (0,2,4):
        for occupied in combinations(range(5),k):
            n3 = sum(i<3 for i in occupied)
            n2 = k-n3
            y = Q(-2*n3+3*n2,6)
            t = Q(int(3 in occupied)-int(4 in occupied),2)
            states[n3,n2,y] += 1
            for i,term in enumerate((t**2,(t+y)**2,y,y**3)):
                traces[i] += term
    assert sorted(states.values()) == [1,1,2,3,3,6]
    assert traces == [2,Q(16,3),0,0]
    assert traces[0]/traces[1] == Q(3,8)
    print("After the constant record is fixed: eigenblock dimensions 3 and 2. Even occupations:")
    for key,dim in sorted(states.items()):
        print(key,"multiplicity",dim)
    print("Weyl trace ratio",traces[0]/traces[1])


if __name__=="__main__":
    main()
