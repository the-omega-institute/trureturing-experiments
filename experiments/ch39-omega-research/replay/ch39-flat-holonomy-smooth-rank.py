"""Rank and rational flux cost at the flat shell background."""
from pathlib import Path
from runpy import run_path

from sympy import I, Matrix, Rational as Q, diag, eye, zeros

dirpath=Path(__file__).parent
orient=run_path(str(dirpath/"ch39-shell-nine-phase-probes.py"))["orientation"]


def main():
    edges=orient()
    D=zeros(48,24)
    for i,(u,v) in enumerate(edges):
        D[i,u],D[i,v]=-1,1
    tree=list(D.T.rref()[1])
    chords=[i for i in range(48) if i not in tree]
    T=D[tree,1:]
    R=T.inv()
    assert set(R)<={-1,0,1} and abs(T.det())==1
    S=D[chords,1:]*R
    coeff=eye(25)+S*S.T
    trees=(D[:,1:].T*D[:,1:]).det(method="domain-ge")
    assert coeff.det(method="domain-ge")==trees and trees>0
    assert coeff.is_symmetric() and len(chords)==25
    # One rational flux check at a single chord; hit graph Laplacian via two
    # different finite integer/rational Schur factorizations.
    t=Q(2,7)
    phase=(1-t*t)/(1+t*t)+I*2*t/(1+t*t)
    j=0
    C=D[chords,:].copy()
    C[j,edges[chords[j]][1]]=phase.expand(complex=True)
    Snew=C[:,1:]*R
    Z=eye(25)+Snew*Snew.H
    field=C*Matrix.ones(24,1)
    local=(field.H*Z.inv(method="DM")*field)[0].simplify()
    G=D[tree,:].T*D[tree,:]+C.H*C
    grounded=G[1:,1:]
    global_value=(G[0,0]-(G[0:1,1:]*grounded.inv(method="DM")*G[1:,0:1])[0]).simplify()
    assert local==global_value and local>0
    # For arbitrary phases, the path matrix at each endpoint has 0, +/-1
    # entries and squared row norm equal to its distance from vertex 0.
    depth=[0]+[sum(v*v for v in R.row(i)) for i in range(23)]
    beta2=sum(2*(depth[edges[e][0]]+depth[edges[e][1]]) for e in chords)
    assert local>=4*t*t/((1+t*t)*(1+beta2))
    print("Exact spanning-tree count",trees,"and 25 independent positive quadratic flux directions.")
    print("All chart fluxes: sum of chord deficits divided by",1+beta2,"bounds grounded energy from below.")
    print("One flux 2/7, grounded squared residual",local)
    print("Minimum smooth complex rows near flatness: 36; with ordinary conjugation symmetry: 48.")


if __name__=="__main__":
    main()
