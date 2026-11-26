"""Two local comparison lengths and the minimal memory of the shell transition."""
from itertools import combinations, permutations
from sympy import Matrix, Rational as Q, eye, ones, symbols, zeros


def main():
    pts = list(permutations((0,1,3,4)))
    pairs = list(combinations(range(4),2))
    split = lambda p:tuple(sorted((p.index(0),p.index(1))))
    def swap(x,y):
        return Matrix([[int(tuple(y if k==x else x if k==y else k for k in p)==q)
                        for q in pts] for p in pts])
    low,high,mid = swap(0,1),swap(3,4),swap(1,3)
    D = Matrix([[26-sum(x*y for x,y in zip(p,q)) for q in pts] for p in pts])
    assert D.applyfunc(lambda x:int(x==1)) == low+high
    assert D.applyfunc(lambda x:int(x==2)) == low*high
    assert D.applyfunc(lambda x:int(x==3)) == zeros(24)
    assert D.applyfunc(lambda x:int(x==4)) == mid
    membership = Matrix([[int(split(p)==s) for s in pairs] for p in pts])
    average = membership.T/4
    E = (eye(24)+low)*(eye(24)+high)/4
    assert E == membership*average and average*membership == eye(6)
    M = average*mid*membership
    assert 4*M == Matrix([[int(len(set(s)&set(t))==1) for t in pairs] for s in pairs])
    distinct = [list(M.row(j)) for j,s in enumerate(pairs) if 0 in s]
    assert {tuple(row) for row in distinct}=={tuple(M.row(i)) for i in range(6)}
    T = Matrix(distinct)
    assert T*T.T == (eye(3)+ones(3))/8
    J = Matrix([[int(set(s).isdisjoint(t)) for t in pairs] for s in pairs])
    assert M*(eye(6)-J)==zeros(6)
    assert (eye(6)-J).rank()==3 and T.rank()==3
    assert M**2*(eye(6)-J)==zeros(6)
    assert average*mid**2*membership==eye(6)
    # Consecutive patch labels identify every microscopic input without a reset.
    histories = set()
    for point in pts:
        first = set(split(point))
        second = {point.index(0),point.index(3)}
        history = (tuple(sorted(first)),tuple(sorted(second)))
        assert history not in histories
        histories.add(history)
        locations = [first & second,first-second,second-first,set(range(4))-(first|second)]
        assert all(len(pos)==1 for pos in locations)
        assert [point[next(iter(pos))] for pos in locations]==[0,1,3,4]
    assert len(histories)==24
    assert average.col_join(average*mid).rank()==11
    energy = 4*(eye(6)-M)
    x = symbols("x")
    assert (energy.charpoly(x).as_expr()-x*(x-4)**3*(x-6)**2).expand()==0
    incidence = Matrix([[int(i in p) for i in range(4)] for p in pairs])
    P = eye(4)-ones(4)/4
    assert T*incidence*P == zeros(3,4)
    assert P*incidence.T*incidence*P == 2*P
    a = symbols("a")
    v = dict(zip((0,1,3,4),(1,a,-a,-1)))
    V = Matrix([[v[k] for k in p] for p in pts])
    assert (average*V*P-(1+a)*incidence*P).applyfunc(lambda b:b.expand())==zeros(6,4)
    # At a threshold of 1/8, complementary output rows still agree to tolerance
    # and all other pairs remain separated for pointwise error strictly below 1/16.
    assert set(max(abs(a-b) for a,b in zip(M.row(i),M.row(j)))
               for i,j in combinations(range(6),2))=={0,Q(1,4)}
    print("Cost 1 gives six disjoint squares; cost 4 first joins them.")
    print("Six equiprobable-patch transitions have three distinguishable output rows, nullity 3 and surjective readout.")


if __name__=="__main__":
    main()
