"""One radial anchor forces a CKM spectator; two endpoints of a shortest arc span all symmetric vertices."""
from itertools import permutations
from sympy import Matrix, Rational as Q, zeros, simplify, sqrt


def main():
    pts = list(permutations((0,1,3,4)))
    a,b = pts[:2]
    assert sum((x-y)**2 for x,y in zip(a,b))==2
    alpha = (sqrt(5)-1)/2
    weight = dict(zip(a,(1,alpha,-alpha,-1)))
    cols = []
    radii = []
    for p in pts:
        z = [weight[t] for t in p[:3]]
        cols.append(Matrix([z[i]*z[j] for i in range(3) for j in range(i,3)]))
        radii.append(tuple(16-sum(s*t for s,t in zip(v,p)) for v in (a,b)))
    def project(values):
        return sum((v*c for v,c in zip(values,cols)),zeros(6,1)).applyfunc(simplify)
    rings = Matrix.hstack(*[project([int(x==d) for x,_ in radii]) for d in sorted({x for x,_ in radii})])
    assert rings.rank()==4
    monomials = [(0,0),(0,2),(1,1),(2,0),(0,4),(1,3)]
    mat = Matrix.hstack(*[project([x**i*y**j for x,y in radii]) for i,j in monomials])
    assert mat.rank()==6
    for i,j in [(1,0),(0,1),(3,0),(2,1),(1,2),(0,3)]:
        assert project([x**i*y**j for x,y in radii])==zeros(6,1)
    print("One anchor radial span 4; two endpoints degrees 0,2,4 give 1,4,6")


if __name__=="__main__":
    main()
