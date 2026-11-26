"""All Wilson masses at the face-action minimum, including the conjugate triplet."""
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path
from flint import arb, acb_mat, ctx


s=run_path(str(Path(__file__).with_name("ch39-shell-gw-symmetric-links.py")))
ctx.prec=240
assert 1-2*arb(3).sqrt()+2<0
# Two cap faces carry 1/4 on all lines. Each line is a spectator in two of six sides.
assert 2*F(1,4)+2*F(5,12)+4*F(-1,12)==1
for sign,interval in ((1,(-.445,-.444,.819,.820)),(-1,(.366,.367,.707,.708)),(0,(-.442,-.441,.532,.533))):
    base=s["model"]["make"](sign,True)
    # At the minimum every link is diagonal. The coordinate cycle relates the three lines.
    graph=s["equivariant_matrix"](base,arb.pi(),arb(0))
    matrix=[[graph[3*i][3*j]-(arb(3)/4 if i==j and i<6 else -arb(3)/4 if i==j else 0)
             for j in range(12)] for i in range(12)]
    polynomial=acb_mat([[(-1 if i<6 else 1)*v for v in row] for i,row in enumerate(matrix)]).charpoly()
    all_roots=polynomial.roots(tol=1e-14)
    real=sorted([z.real for z in all_roots if z.imag.contains(0)],key=lambda z: float(z.mid()))
    assert len(real)==2 and interval[0]<real[0]<interval[1] and interval[2]<real[1]<interval[3]
    def effective_index(m):
        current=[row.copy() for row in matrix]
        for i in range(12):
            current[i][i]+=m if i<6 else -m
        return 6-s["model"]["count"](current,arb(0))
    assert [effective_index(m) for m in (-2,(real[0]+real[1]).mid()/2,2)]==[0,1,0]
    print("center sign",sign,"index one between",*[v.str(13) for v in real],"; zero outside")
