"""曾登记的两枚读弱零点在交换或共用时的相邻奇异值比。"""
from pathlib import Path
from runpy import run_path
from flint import acb, acb_mat, arb, ctx
R=run_path(str(Path(__file__).with_name("ch39-moving-preparation-quark-shape.py")))
lib=R["lib"]
g=lib("three-ray-chiral-boundary")
adj=lib("native-light-projector-pump")["adjoint"]
spectrum=lib("adjacent-pvm-transport")["spectrum"]
sv=lib("first-zero-relative-quark-spectrum")["singular_values"]
ctx.prec=320
W=R["orthonormal_sources"](g,adj)
port=[g["soft"]["pole"],*g["boundary_sites"]()]
table=lib("native-rh-ray-pvm")["gram_cholesky"](R["ROOT"])
G=acb_mat([[table[j][port.index(i)] if i in port else 0 for j in (1,2,3)] for i in range(19)])
print("ports",port,"s(G)",*[z.str(15) for z in sv(G)])
M={}
raw={}
for q in (-4,2):
    for n in (1,2086602,2040244):
        e,U=spectrum(g,q,n)
        V=acb_mat([[U[i,j] for j in range(3)] for i in range(19)])
        filt=acb_mat([[1/(1+(e[i]/9)**512) if i==j else 0 for j in range(3)] for i in range(3)])
        M[q,n]=filt*adj(V)*W
        raw[q,n]=filt*adj(V)*G
        print(q,n,"s(A)",*[z.str(12) for z in sv(M[q,n])])
for n,m in ((1,1),(2086602,2086602),(2040244,2040244),(2086602,2040244),(2040244,2086602)):
    A,B=M[-4,n],M[2,m]
    s=sv(A*B.inv())
    t=sv(B*A.inv())
    h=sv(raw[-4,n]*raw[2,m].inv())
    assert all(x.overlaps(1/y) for x,y in zip(s,reversed(t)))
    assert all(x.overlaps(y) for x,y in zip(s,h))
    print("heights",n,m,"up/down adjacent", (s[1]/s[0]).str(15),(s[2]/s[1]).str(15))
