"""Five shell angle ratios determine a real three-axis response up to scale."""
from itertools import combinations, product
from sympy import Matrix, Rational


def form(p,q):
    return [p[i]*q[i] for i in range(3)]+[p[i]*q[j]+p[j]*q[i] for i,j in ((0,1),(0,2),(1,2))]


def main():
    points=[p for p in product(range(5),repeat=3) if sum(i*i for i in p)==17]
    for numerator,goal in [(33,200),(28,227)]:
        rows=[]
        for p,q in combinations(points,2):
            dot=sum(i*j for i,j in zip(p,q))
            if 17+dot==numerator:
                rows.append([34*x-dot*(y+z) for x,y,z in zip(form(p,q),form(p,p),form(q,q))])
        assert len(rows)==6
        for chosen in combinations(rows,5):
            m=Matrix([*chosen,[1,1,1,0,0,0]])
            assert m.det()!=0 and m*Matrix([1,1,1,0,0,0])==Matrix([0,0,0,0,0,3])
        inv=Matrix(rows[:5]+[[1,1,1,0,0,0]]).inv()
        v=[1156*sum(abs(entry) for entry in row[:5]) for row in inv.tolist()]
        bound=max(v[0]+v[3]+v[4],v[1]+v[3]+v[5],v[2]+v[4]+v[5])
        assert bound<goal
        print("Any five out of six arcs",numerator,"/34 fix isotropy. First five conditioning factor",bound)
    # The absolute scale is supplied by the trace or one sum-and-difference total.
    # In the latter case |trace(H)/3 - measured scale| <= ||H-trace(H)I/3||.
    eps=Rational(1,2**112)
    for factor in (200,227):
        assert (2*factor*eps*(1+eps)/(1-2*factor*eps)+eps)<Rational(1,2**102)
    print("112-bit ratios and absolute channel sum give 102-bit matrix precision.")


if __name__=="__main__":
    main()
