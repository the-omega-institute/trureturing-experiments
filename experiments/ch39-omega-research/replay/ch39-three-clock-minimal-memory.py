"""Minimal Hermitian memory at the three clock roots, in rational weighted coordinates."""
from pathlib import Path
from runpy import run_path

from flint import fmpq, fmpq_mat

REFINE = run_path(str(Path(__file__).with_name("ch39-refined-mass-obstruction.py")))


def identity(n):
    return fmpq_mat([[int(i==j) for j in range(n)] for i in range(n)])


def coefficients():
    nodes, edges = REFINE["geometry"]()
    m, n = 19, len(nodes)-19
    injection, internal = fmpq_mat(n,m), fmpq_mat(n,n)
    for i,j,w,charge in edges:
        if j < m:
            continue
        assert charge==0
        weight = fmpq(w.numerator,w.denominator)
        if i < m:
            injection[j-m,i] = weight
        else:
            internal[i-m,j-m] = internal[j-m,i-m] = weight
    advanced = internal*injection
    columns = fmpq_mat([row+other for row,other in zip(injection.tolist(),advanced.tolist())])
    assert injection.rank()==19 and columns.rank()==38
    r = (identity(n)-internal**3).inv()
    g0 = injection.transpose()*r*injection
    g1 = injection.transpose()*r*advanced
    g2 = advanced.transpose()*r*advanced
    return internal, columns, r, g0, g1, g2


def verify():
    D, X, r, g0, g1, g2 = coefficients()
    H = X.transpose()*r*X
    K = H.solve(X.transpose()*r*D*X)
    E = fmpq_mat([[int(i==j) for j in range(19)] for i in range(38)])
    second = fmpq_mat([[int(i==j+19) for j in range(19)] for i in range(38)])
    assert K*E==second and H*K==K.transpose()*H
    assert [E.transpose()*H*K**j*E for j in range(3)]==[g0,g1,g2]
    gram = H*(identity(38)-K**3)
    assert gram==gram.transpose() and gram*K==K.transpose()*gram
    k = g0.solve(g1)
    shift = g2-g1*k
    assert shift==shift.transpose() and shift.rank()==19
    small_gram = g0*(identity(19)-k**3)
    assert small_gram==small_gram.transpose() and small_gram*k==k.transpose()*small_gram
    assert [g0, g0*k, shift+g0*k*k]==[g0,g1,g2]
    print("Old block fixed: auxiliary dimension",H.nrows(),
          "old block adjusted: dimension",g0.nrows(),"shift rank",shift.rank())
    print("Trace of required old-block shift",float(sum(shift[i,i] for i in range(19))))


if __name__ == "__main__":
    verify()
