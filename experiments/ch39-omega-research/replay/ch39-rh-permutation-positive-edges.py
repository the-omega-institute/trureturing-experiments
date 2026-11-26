"""Signed-edge cost of the two-arc charged metric, including hidden-node reductions."""
from itertools import permutations
from pathlib import Path
from runpy import run_path

from flint import arb,ctx
from sympy import Matrix, Rational as Q, eye, ones, symbols, zeros

arc = run_path(str(Path(__file__).with_name("ch39-rh-internal-arc-calibration.py")))


def class_sums():
    states=list(permutations(range(4)))
    mats={key:zeros(24) for key in ("T","D","C","F")}
    for g in states:
        moved=sum(i!=u for i,u in enumerate(g))
        if not moved:
            continue
        parity=sum(g[i]>g[j] for i in range(4) for j in range(i+1,4))%2
        label={2:"T",3:"C",4:"F" if parity else "D"}[moved]
        for i,s in enumerate(states):
            mats[label][i,states.index(tuple(s[g[j]] for j in range(4)))]+=1
    return [mats[key]/sum(mats[key][0,:]) for key in ("T","D","C","F")]


def main():
    ctx.prec, ctx.cap=1536,3
    data=arc["fixed"]["inputs"]()
    d=arc["argument_for"](arb(14)/17,data)/(3*arc["argument_for"](arb(33)/34,data))
    assert d>3
    trans, double, triple, four=class_sums()
    ident=eye(24)
    assert four**2==ident/6+double/6+Q(2,3)*triple
    # Columns: trivial, sign, ordinary standard, sign-standard, pair representation.
    lap=Matrix([[0,2,Q(2,3),Q(4,3),1],[0,0,Q(4,3),Q(4,3),0],
                [0,0,1,1,Q(3,2)],[0,2,Q(4,3),Q(2,3),1]])
    assert lap[:,1]-3*lap[:,3]==Matrix([-2,-4,-3,0])
    a,D=symbols("a D")
    coeff=Matrix([0,-(D-3+2*a)/4,0,(D-a)/2])
    assert (a*ones(1,5)+coeff.T*lap-Matrix([[a,D,(D+3-a)/3,1,(a+D)/2]])).applyfunc(lambda e:e.expand())==zeros(1,5)
    # Verify the full 24-dimensional operator identity for these five characters.
    sizes=(1,1,3,3,2)
    traces=Matrix([[1]*5]).col_join(ones(4,5)-lap)
    parts=[]
    for j,dimension in enumerate(sizes):
        projector=dimension**2*sum((size*traces[i,j]*m for i,(size,m) in
                                  enumerate(zip((1,6,3,8,6),(ident,trans,double,triple,four)))),zeros(24))/24
        assert projector**2==projector and projector.trace()==dimension**2
        assert all((ident-m)*projector==lap[i,j]*projector for i,m in enumerate((trans,double,triple,four)))
        parts.append(projector)
    assert sum(parts,zeros(24))==ident
    print("At constant and colored mean both one, minimum positive off-diagonal row sum",((d-1)/4).str(15))
    print("If constant is merely positive, sharp infimum",((d-3)/4).str(15))
    print("Other regular eigenvalues (ordinary standard, pair representation)",((d+2)/3).str(14),((d+1)/2).str(14))


if __name__=="__main__":
    main()
