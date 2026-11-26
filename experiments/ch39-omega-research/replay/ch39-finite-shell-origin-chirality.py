"""A complex of four-point orbit sums and their local compatibility equations."""
from itertools import permutations

from sympy import Matrix, eye, symbols, zeros


def compose(p,q):
    return tuple(p[i] for i in q)


def cosets(group,subgroup):
    return sorted({tuple(sorted(compose(p,h) for h in subgroup)) for p in group})


def main():
    group = list(permutations(range(4)))
    unit = group[0]
    klein = {p for p in group if p==unit or (all(p[i]!=i for i in unit) and compose(p,p)==unit)}
    cyclic = set()
    for p in group:
        p2 = compose(p,p)
        if p2!=unit and compose(p2,p2)==unit:
            cyclic.add(tuple(sorted((unit,p,p2,compose(p2,p)))))
    subgroups = [tuple(sorted(klein))]+sorted(cyclic)
    assert len(subgroups)==4 and all({h[0] for h in sub}==set(unit) for sub in subgroups)
    patches = [patch for sub in subgroups for patch in cosets(group,sub)]
    assert len(patches)==len(set(patches))==24
    B = Matrix([[int(p in patch) for p in group] for patch in patches])
    rows = []
    for i,sub in enumerate(subgroups[1:],1):
        joint = {compose(k,h) for k in klein for h in sub}
        assert len(joint)==8 and all(compose(x,y) in joint for x in joint for y in joint)
        for big in cosets(group,joint):
            row = [0]*24
            for j,patch in enumerate(patches):
                if j//6 in (0,i) and set(patch)<=set(big):
                    row[j] = -1 if j//6==0 else 1
            assert row.count(-1)==row.count(1)==2
            rows.append(row)
    A = Matrix(rows)
    assert A*B==zeros(9,24)
    H = Matrix([[1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]])/2
    X = Matrix.vstack(*[(H.T*H.extract(p,range(3))).reshape(1,9) for p in group])
    assert B*X==zeros(24,9) and X.T*X==8*eye(9)
    assert B.rank()==15 and A.rank()==9
    # Even cells are the 24 microscopic records and 9 consistency records.
    Q = B.row_join(A.T)
    x = symbols("x")
    assert (B.T*B).charpoly(x).as_expr().factor()==x**9*(x-4)**10*(x-10)**4*(x-16)
    assert (A*A.T).charpoly(x).as_expr().factor()==(x-8)*(x-5)**4*(x-2)**4
    assert (Q*Q.T).charpoly(x).as_expr().factor()==(x-16)*(x-10)**4*(x-8)*(x-5)**4*(x-4)**10*(x-2)**4
    print("24 -> 24 -> 9, ranks 15 and 9, even kernel 3 x 3, odd kernel zero, squared gap 2.")


if __name__=="__main__":
    main()
