"""Signed position-label marginals and the independent three-outcome weak register."""
from itertools import combinations, permutations

from sympy import Matrix, diag, eye, ones, zeros


def main():
    pts = list(permutations(range(4)))
    parity = Matrix([(-1)**sum(p[i]>p[j] for i,j in combinations(range(4),2)) for p in pts])
    H = Matrix([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]])/2
    X = Matrix.vstack(*[(H.T*H.extract(p,range(3))).reshape(1,9) for p in pts])
    signed = diag(*parity)*X
    counts = Matrix([[parity[j]*int(p[i]==k) for j,p in enumerate(pts)]
                     for i in range(4) for k in range(4)])
    assert counts.rank()==10
    assert counts.col_join(signed.T).col_join(parity.T).rank()==10
    B = signed.T.col_join(ones(1,24)).col_join(parity.T)
    Q = diag(B,ones(1,3))
    assert Q*Q.T == diag(8*eye(9),24,24,3)
    assert B*X==zeros(11,9)
    # The unsigned pairing-matrix coefficients give four further independent
    # microscopic directions. The nonconstant output register supplies two.
    pairs = [(0,1),(0,2),(0,3)]
    cuts = [frozenset((frozenset(p),frozenset(set(range(4))-set(p)))) for p in pairs]
    E = Matrix([[1,0],[0,1],[-1,-1]])
    weak = []
    for p in pts:
        dest = [frozenset(frozenset(p[i] for i in side) for side in cut) for cut in cuts]
        R = Matrix([[int(r==s) for s in cuts] for r in dest])
        weak.append((E.T*R*E).reshape(1,4))
    W = Matrix.vstack(*weak)
    kernel = diag(X.row_join(W),E)
    assert kernel.rank()==15 and Q*kernel==zeros(12,15)
    print("Parity-marginal model: 27 -> 12, 15 even zero modes, no odd, squared gap 3.")
    print("Nonzero squared singular values: 8 times nine, 24 times two, 3 times one.")


if __name__=="__main__":
    main()
