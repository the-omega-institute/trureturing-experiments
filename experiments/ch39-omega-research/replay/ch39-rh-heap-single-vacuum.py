"""Exact Hessians and charged-mass compressions in each four-cycle vacuum."""
from itertools import permutations
from pathlib import Path
from runpy import run_path

from sympy import Matrix, Rational, diag, eye, zeros

zero=run_path(str(Path(__file__).with_name("ch39-rh-heap-zero-fields.py")))
states,_,four,cosets,cid,even,odd=zero["quotient"]()
compose,parity=zero["compose"],zero["parity"]


def main():
    h=((1,1,1),(1,-1,-1),(-1,1,-1),(-1,-1,1))
    cols=Matrix([[parity(p)*sum(h[r][i]*h[p[r]][j] for r in range(4))//4
                  for i in range(3) for j in range(3)]+[1,parity(p)] for p in states])
    assert cols.T*cols==diag(8*eye(9),24,24)
    index={p:i for i,p in enumerate(states)}
    base_spectrum={0:1,Rational(8,3):9,8:13,16:1}
    triad_spectrum={0:2,Rational(8,3):13,Rational(16,3):3,8:5,Rational(32,3):1}
    for match in [None,*permutations(odd)]:
        lines=[0]*6
        if match:
            for i,(v,w) in enumerate(zip(even,match)):
                lines[v]=lines[w]=i
        for axis in (1,2,3):
            rows=[]
            for p in states:
                for g in four:
                    q=compose(p,g)
                    r=compose(q,g)
                    a,b=lines[cid[p]],lines[cid[q]]
                    orientation=1 if a==b or (a^b)==axis else -1
                    row=[0]*24
                    row[index[p]],row[index[q]],row[index[r]]=orientation,-1-orientation,1
                    rows.append(row)
            gradient=Matrix(rows)
            curvature=gradient.T*gradient/6
            singlet=Rational(80,9) if match else 16
            expected=diag(Rational(64,3)*eye(9),0,24*singlet)
            assert cols.T*curvature*cols==expected
            assert curvature.eigenvals()==(triad_spectrum if match else base_spectrum)
    print("All 21 pattern-axis calculations: colored 8/3 with multiplicity nine, no singlet mixing.")
    print("Full tangent kernel dimensions: three for the collinear fields, six for the triad fields.")
    print("Positive Hessian gap in both types: 8/3.")


if __name__=="__main__":
    main()
