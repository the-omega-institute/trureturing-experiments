"""An oriented-arc flavor projector commuting with the finite Wilson kernel and proton protection charges."""
from itertools import permutations
from pathlib import Path
from runpy import run_path

from sympy import Matrix, Rational as Q, eye, sqrt, zeros, simplify

shell = run_path(str(Path(__file__).with_name("ch39-finite-shell-first-triplet.py")))


def main():
    pts = list(permutations((0,1,3,4)))
    L = shell["matrix"]({1:1,2:1,4:1})
    def exchange(i,j):
        def to(p):
            v = list(p)
            v[i],v[j] = v[j],v[i]
            return tuple(v)
        return Matrix([[int(to(p)==q) for q in pts] for p in pts])
    P = (eye(24)+exchange(0,1))*(eye(24)+exchange(2,3))/4
    assert P*P==P and P.T==P and P*L==L*P and P.rank()==6
    alpha = (sqrt(5)-1)/2
    v = dict(zip(pts[0],(1,alpha,-alpha,-1)))
    V = Matrix([[v[p[0]]+v[p[i]] for i in (1,2,3)] for p in pts])
    assert (P*V-V*Matrix.diag(1,0,0)).applyfunc(simplify)==zeros(24,3)
    # Columns: color dimension, weak dimension, signed SM hypercharge, X.
    rows = []
    for ell in (-9,0,0):
        rows.extend([(3,2,Q(1,6),1),(3,1,-Q(2,3),-1),(3,1,Q(1,3),-1),
                     (1,2,-Q(1,2),ell),(1,1,1,-ell),(1,1,0,-ell)])
    assert sum(w*x for c,w,y,x in rows if c==3)==0
    assert sum(c*x for c,w,y,x in rows if w==2)==0
    for h,k in ((0,1),(0,3),(1,2),(2,1)):
        assert sum(c*w*y**h*x**k for c,w,y,x in rows)==0
    # The unnormalized one-loop hypercharge/X kinetic mixing trace is nonzero.
    assert sum(c*w*y*x for c,w,y,x in rows)==24
    print("Arc projector rank 6 on 24 sites, rank 1 on the chiral triplet; X anomalies vanish, Tr(YX)=24.")


if __name__=="__main__":
    main()
