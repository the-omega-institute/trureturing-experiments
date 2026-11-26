"""Produce the four polarity channels on the original finite permutation shell."""
from itertools import combinations,permutations,product
from pathlib import Path
from runpy import run_path

from sympy import Matrix, Rational, diag, eye, zeros

directory=Path(__file__).parent
zero=run_path(str(directory/"ch39-rh-heap-zero-fields.py"))
P,V,four,cosets,cid,even,odd=zero["quotient"]()
compose=zero["compose"]
cls=run_path(str(directory/"ch39-rh-permutation-positive-edges.py"))


def main():
    index={p:i for i,p in enumerate(P)}
    translations={}
    for g in V:
        r=zeros(24)
        for i,p in enumerate(P):
            r[i,index[compose(p,g)]]=1
        translations[g]=r
    chars=[[1,1,1,1]]+[[1]+[1 if i==j else -1 for i in range(3)] for j in range(3)]
    proj=[sum((s*translations[g] for s,g in zip(v,V)),zeros(24))/4 for v in chars]
    for i,j in product(range(4),repeat=2):
        assert proj[i]*proj[j]==(proj[i] if i==j else zeros(24))
        signed=sum((4*int(g==V[0])-1)*chars[i][k]*chars[j][k] for k,g in enumerate(V))
        assert signed==4*int(i!=j)
    assert sum(proj,zeros(24))==eye(24) and all(p.rank()==6 for p in proj)
    _,double,_,four_mean=cls["class_sums"]()
    K=6*eye(24)-8*four_mean+2*double
    assert all(K*m==m*K for m in proj)
    for i,m in enumerate(proj):
        basis=Matrix.hstack(*m.columnspace())
        compressed=(basis.T*basis).inv()*basis.T*K*basis
        expected={0:1,8:4,16:1} if i==0 else {Rational(8,3):3,8:3}
        assert compressed.eigenvals()==expected
    # The six ordered character labellings are exactly the conjugation orbit.
    inverses={p:tuple(p.index(i) for i in range(4)) for p in P}
    images={tuple(V.index(compose(compose(p,g),inverses[p])) for g in V[1:]) for p in P}
    assert images==set(permutations((1,2,3)))
    print("Four coset characters realized by relative comparisons, six copies of each in 24 dimensions.")
    print("Same-versus-different relative label recovers all the raw phase cross terms.")
    print("Trivial sector spectrum 0,8,16 (1,4,1); each nontrivial sector 8/3 and 8 (3,3).")


if __name__=="__main__":
    main()
