"""Checks for La Piana, arXiv:2512.08797v3, Section 4.1 (Question 1.3 for C_5): the internal-dimension-2 block matrix
B_ij = (1/(5n))[[5+a_(i-j), b_(i+j)],[b_(i+j), 5-a_(i-j)]] is a graph quantum magic square of C_n (PSD blocks, magic relations, block
graph relation) for odd n = 5..13, and the rank-one data used by the separation argument hold."""
import sympy as sp
ok=True
for n in (5,7,9,11,13):
    a=lambda r: {0:3,1:-3}.get(r%n,0); b=lambda s: {0:4,1:-4}.get(s%n,0)
    B=[[sp.Matrix([[5+a(i-j), b(i+j)],[b(i+j), 5-a(i-j)]])/(5*n) for j in range(n)] for i in range(n)]
    psd=all(B[i][j][0,0]>=0 and B[i][j][1,1]>=0 and B[i][j].det()>=0 for i in range(n) for j in range(n))
    rows=all(sum((B[i][j] for j in range(n)),sp.zeros(2))==sp.eye(2) for i in range(n))
    cols=all(sum((B[i][j] for i in range(n)),sp.zeros(2))==sp.eye(2) for j in range(n))
    adj=lambda k,l: (k-l)%n in (1,n-1)
    graph=all(sum((B[i][k] for k in range(n) if adj(k,j)),sp.zeros(2))==sum((B[k][j] for k in range(n) if adj(k,i)),sp.zeros(2)) for i in range(n) for j in range(n))
    F=lambda r,s: sp.Matrix([[5+a(r),b(s)],[b(s),5-a(r)]])/(5*n)
    # rank-one facts used by the separation: F00,F01 rank one with distinct ranges; F10,F11 likewise; F00 != F10 top-left
    r1=all(F(r,s).rank()==1 for r in (0,1) for s in (0,1))
    distinct=F(0,0).columnspace()[0].cross(sp.Matrix([0,0,0])) is None if False else True
    sep=(sp.Matrix.hstack(F(0,0),F(0,1)).rank()==2) and (sp.Matrix.hstack(F(1,0),F(1,1)).rank()==2) and F(0,0)[0,0]!=F(1,0)[0,0]
    good=psd and rows and cols and graph and r1 and sep
    print(n,'PSD',psd,'rows',rows,'cols',cols,'graph',graph,'rank-one',r1,'separation data',sep)
    ok&=good
print('RESULT','PASS' if ok else 'FAIL')
