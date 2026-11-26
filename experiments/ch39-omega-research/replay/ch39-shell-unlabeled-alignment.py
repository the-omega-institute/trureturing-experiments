"""Certify the optimal injective alignment over all rotations of the two unlabeled shells."""
from itertools import product
from pathlib import Path
from flint import fmpz_mat,arb
import runpy

geom=runpy.run_path(str(Path(__file__).with_name('ch39-shell-family-alignment.py')))
source,target=geom["sites"],geom["shell"](867)
order=(0,5,2,3,1,4,6,7,8)
S=[source[i] for i in order]
magnitude=geom["den"]
# Trace threshold for candidate preservation. Below the optimum to isolate nearby alternatives.
cut=1085
def score(x,y):
    c=[[a*b for b in y] for a in x]
    tr=sum(c[i][i] for i in range(3))
    z=[c[(i+1)%3][(i+2)%3]-c[(i+2)%3][(i+1)%3] for i in range(3)]
    return fmpz_mat([[tr,*z]]+[[z[i],*[c[i][j]+c[j][i]-int(i==j)*tr for j in range(3)]] for i in range(3)])
components=[[score(x,y) for y in target] for x in S]
def chord(x,y,r):
    return (sum((u-v)**2 for u,v in zip(x,y))/arb(r)).sqrt()
D=[[chord(x,y,17) for y in S] for x in S]
E=[[chord(x,y,867) for y in target] for x in target]
slack=(2*(18-2*cut/magnitude)).sqrt()
# If the total squared point error is <= eps, every chord discrepancy <= sqrt(2 eps).
pairs={(i,j,k):{l for l in range(19) if abs(D[i][j]-E[k][l])<slack}
       for i in range(9) for j in range(i) for k in range(19)}
assert all(abs(abs(d-e)-slack)>0 for row in D for d in row for line in E for e in line)
unit=fmpz_mat([[int(i==j) for j in range(4)] for i in range(4)])
def pd(A):
    # The first three leading principal minors have cheaper exact formulas.
    a=A[0,0]
    if a<=0:
        return False
    b=a*A[1,1]-A[1,0]**2
    if b<=0:
        return False
    c=b*A[2,2]-a*A[2,1]**2+2*A[1,0]*A[2,0]*A[2,1]-A[1,1]*A[2,0]**2
    return c>0 and A.det()>0
assert magnitude<arb(121405)/1000
candidates=[]
nodes=[0]*10
def visit(seq,K):
    depth=len(seq)
    nodes[depth]+=1
    threshold=cut*1000-(9-depth)*121405
    if pd(threshold*unit-1000*K):
        return
    if depth==9:
        candidates.append((seq,K))
        return
    options=set(range(19))-set(seq)
    for j,l in enumerate(seq):
        options.intersection_update(pairs[depth,j,l])
    for l in sorted(options):
        visit(seq+(l,),K+components[depth][l])
visit((),fmpz_mat(4,4))
by_poly={}
for seq,K in candidates:
    by_poly.setdefault(tuple(map(int,K.charpoly().coeffs())),[]).append(seq)
print("full assignments at or above cutoff",cut,len(candidates),"surviving polynomials",len(by_poly))

from flint import fmpz_poly
ranked=[]
for poly,seqs in by_poly.items():
    roots=fmpz_poly(list(poly)).complex_roots()
    assert all(root.imag.contains(0) for root,n in roots)
    ordered=sorted((root.real for root,n in roots),key=lambda z:float(z.mid()))
    assert len(ordered)==1 or ordered[-1]>ordered[-2]
    ranked.append((ordered[-1],seqs,poly))
ranked.sort(key=lambda item:float(item[0].mid()),reverse=True)
assert all(x[0]>y[0] for x,y in zip(ranked,ranked[1:]))
winners=ranked[0][1]
expected=set()
for shift in range(3):
    expected.add(tuple(target.index(tuple(geom["digits"][point[(i+shift)%3]] for i in range(3)))
                       for point in S))
assert expected==set(winners)
runner=ranked[1][0]
gap=2*(1089-runner)/magnitude
assert len(candidates)==84 and len(ranked)==7 and len(ranked[1][1])==9
assert arb(".02227")<gap<arb(".02229")
print("number of best / second injections",len(winners),len(ranked[1][1]))
print("second-best correlation polynomial low-to-high",ranked[1][2])
print("minimum sum", (18-2*1089/magnitude).str(15),"gap",gap.str(15))
