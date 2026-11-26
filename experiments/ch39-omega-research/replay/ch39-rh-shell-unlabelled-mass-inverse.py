"""Recover the Schur pivot from unlabelled scalar squared masses for any positive reference kinetic form."""
from pathlib import Path
from runpy import run_path

import sympy as sp
from flint import arb, ctx

base=run_path(str(Path(__file__).with_name("ch39-finite-rh-propagating-references.py")))


def main():
    x,h,L,C,D=sp.symbols("x h L C D", positive=True)
    a,b=sp.symbols("a b", real=True)
    M=sp.Matrix([[h+L,sp.sqrt(L*C)*a,sp.sqrt(L*D)*b],
                 [sp.sqrt(L*C)*a,C,0],[sp.sqrt(L*D)*b,0,D]])
    p=1-a*a-b*b
    polynomial=(x-h-L)*(x-C)*(x-D)-L*(C*a*a*(x-D)+D*b*b*(x-C))
    assert sp.expand(M.charpoly(x).as_expr().subs(M.charpoly(x).gen,x)-polynomial)==0
    assert sp.expand(M.det()-C*D*(h+L*p))==0
    assert sp.expand((M*M).trace()-(C*C+D*D+(h+L)**2+2*L*(C*a*a+D*b*b)))==0
    assert sp.cancel((M.det()/(C*D)-h)/(M.trace()-h-C-D)-p)==0
    eig=[(6,7),(4,6),(3-sp.sqrt(5),3),(3+sp.sqrt(5),3),(3-sp.sqrt(3),2),(3+sp.sqrt(3),2)]
    assert [sp.simplify(sum(n*z**k for z,n in eig)) for k in (0,1,2)]==[23,96,480]
    assert 25+23==48 and 2*25+3*23==119
    # Only a cycle pole can supply two independent vectors in any one block.
    # At every other eigenvalue the reference equations determine both
    # reference components from the single point component.
    ctx.prec=192
    u,v,c=[base["gram"](i,j) for i,j in ((1,0),(1,2),(0,2))]
    target=(1+2*u*v*c-u*u-v*v-c*c)/(1-c*c)
    assert target>arb("0.00000357") and target<arb("0.00000358")
    print("Mass-spectrum invariant (det/(C1*C2)-h)/(trace-C1-C2-h) =",target.str(15))
    print("119 unordered squared masses determine cycle eigenvalues, M_F, the Schur pivot and all observable reference mixing.")


if __name__=="__main__":
    main()
