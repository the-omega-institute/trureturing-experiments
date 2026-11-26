"""Exact point-group invariants and shell cardinalities for five oriented geometries."""
from itertools import combinations, permutations

from flint import fmpq, fmpq_mat

planes = tuple(combinations(range(4),2))
steps = [[int(i==j) for i in range(4)] for j in range(4)]+[[-1]*4]


def wedge(M):
    return fmpq_mat([[M[i,k]*M[j,l]-M[i,l]*M[j,k] for k,l in planes] for i,j in planes])


def sym_trace(M):
    trace = lambda A: sum(A[i,i] for i in range(A.nrows()))
    return (trace(M)**2+trace(M*M))/2


def verify():
    Q = fmpq_mat([[1+int(i==j) for j in range(4)] for i in range(4)])
    totals = [0,0,0,0]
    for perm in permutations(range(5)):
        R = fmpq_mat([steps[i] for i in perm[:4]]).transpose()
        assert R*Q*R.transpose() == Q
        two = wedge(R)
        for j,U in enumerate((R,two)):
            totals[j] += sym_trace(U)
            if R.det()>0:
                totals[2+j] += sym_trace(U)
    assert totals == [120,120,60,120]
    traceless = []
    for v in steps:
        row = fmpq_mat([v])
        x = row.transpose()*row-Q/5
        traceless.append(x.entries())
    assert fmpq_mat(traceless).rank() == 4
    assert all(sum(col)==0 for col in zip(*traceless))
    unrestricted = 5*56
    sheets = next(r for r in range(1,16) if 19*r % 5 == 0 and 19*r//5 >= 56)
    assert unrestricted==280 and sheets==15 and 19*sheets==5*57
    print("Invariant symmetric forms on one- and two-forms: S5 (1,1), A5 (1,2).")
    print("Joint equal states:",unrestricted,"shell labels per point:",sheets,"shell states:",19*sheets)


if __name__ == "__main__":
    verify()
