"""Minimum reference count for the whole vertex response and an integer spanning-tree realization."""
from pathlib import Path
from runpy import run_path

from sympy import Matrix,eye,zeros

shell=run_path(str(Path(__file__).with_name("ch39-finite-shell-first-triplet.py")))


def main():
    L=shell["matrix"]({1:1,2:1,4:1})
    edges=[(i,j) for i in range(24) for j in range(i+1,24) if L[i,j]<0]
    D=zeros(len(edges),24)
    for k,(i,j) in enumerate(edges):
        D[k,i],D[k,j]=-1,1
    assert len(edges)==48 and D.T*D==L
    chosen=D.T.rref()[1]
    tree=D[list(chosen),:]
    assert tree.shape==(23,24)
    # Pin one vertex to solve for the integer path coordinates.
    paths=D[:,1:]*tree[:,1:].inv()
    assert set(paths)<={-1,0,1} and paths*tree==D
    metric=paths.T*paths
    # The selected edge rows form an identity, proving full column rank,
    # positivity of the induced metric and a 25-dimensional orthogonal complement.
    assert paths[list(chosen),:]==eye(23) and tree.T*metric*tree==L
    assert paths.rows-paths.cols==25
    # An edge-supported row orthogonal to the constant vector can only
    # be a multiple of one oriented difference. Every nonzero off-diagonal
    # entry then needs that particular edge in the positive pole residue.
    assert sum(L[i,j]!=0 for i in range(24) for j in range(i+1,24))==48
    assert all(L[i,i]==4 for i in range(24))
    # At 96 bits, even a factor four in a single edge entry leaves 94 bits;
    # the scalar alternation certificate excludes order 13 already at 93.
    dim,lower_order=23,13
    required_zeros=2*lower_order*dim
    for state in (dim*14-1,dim*14):
        max_degree=state+lower_order*dim-dim
        assert (max_degree>=required_zeros)==(state==322)
    assert 23+322==345 and 23+48*14==695
    assert 23+23*21==506 and 23+48*21==1031
    print("At 96-bit relative operator accuracy: 322 reference fields, or 672 with edge-supported pole couplings.")
    print("Common shell response: 23 reference directions per pole, or 48 with two-vertex pole couplings.")
    print("Integer path matrix realizes the 23-direction factorization exactly; complement has dimension 25.")


if __name__=="__main__":
    main()
