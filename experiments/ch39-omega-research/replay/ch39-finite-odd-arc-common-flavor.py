"""A binary shortest-arc reference gives complementary flavor readings of one scalar profile."""
from itertools import permutations
from sympy import Matrix, Rational as Q, diag, eye, simplify, sqrt, zeros


def main():
    pts = list(permutations((0,1,3,4)))
    a,b = pts[:2]
    alpha = (sqrt(5)-1)/2
    lookup = dict(zip(a,(1,alpha,-alpha,-1)))
    cols,signs = [],[]
    for p in pts:
        z = [lookup[x] for x in p]
        # Real Hadamard axes for the sum-zero subspace, with common squared
        # zero-mode norm 8 per unnormalized unit feature in the full 24x3 frame.
        v = Matrix([z[0]+z[1],z[0]+z[2],z[0]+z[3]])
        P = (v*v.T/(8*(2+2*alpha**2))).applyfunc(simplify)
        cols.append(Matrix([P[i,i] for i in range(3)]+[sqrt(2)*P[i,j] for i,j in ((0,1),(0,2),(1,2))]))
        value = sum(x*(i-j) for x,i,j in zip(p,a,b))
        assert value!=0
        signs.append(1 if value>0 else -1)
    F = Matrix.hstack(*cols)
    gram = (F*F.T).applyfunc(simplify)
    off = diag(*([Q(17,160)]*3+[Q(1,80)]*3))
    one = Matrix([1,1,1,0,0,0])
    assert gram == off+(Q(1,8)-Q(17,160))*one*one.T/3
    assert (F*diag(*signs)*F.T).applyfunc(simplify)==zeros(6)
    up = [Q(s) for s in (".00000615",".00311",".8616")]
    u = sum(t*t for t in up)
    e = sum(up)**2/3
    down = sum(Q(s)**2 for s in (".0000135",".000268",".01401"))
    n2 = 8*e+Q(160,17)*(u-e)+80*down
    assert n2 < Q("2.579")**2
    # Color-singlet s-channel in the high-energy approximation with frozen columns.
    assert 3*n2/8 < 16*Q(314,100)*Q(".05")
    print("Two six-dimensional row spaces are orthogonal; common-profile norm below 2.579.")


if __name__=="__main__":
    main()
