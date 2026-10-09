import itertools, sympy as sp
from collections import Counter
def admissible(w):
    n=len(w)
    if n<=1: return all(c==0 for c in w)   # M_0={λ}, M_1={0}
    if all(c==1 for c in w): return False
    if all(c==0 for c in w): return True
    # rotate so that w starts at the beginning of a 1-run (previous char 0)
    k=next(i for i in range(n) if w[i]==1 and w[i-1]==0)
    v=w[k:]+w[:k]
    # parse runs
    runs=[]; i=0
    while i<n:
        j=i
        while j<n and v[j]==v[i]: j+=1
        runs.append((v[i],j-i)); i=j
    # runs alternate starting with 1, ending with 0
    for a in range(0,len(runs),2):
        if runs[a+1][1] <= runs[a][1]: return False
    return True
def N(n):
    words=[w for w in itertools.product([0,1],repeat=n) if admissible(list(w))]
    S=set(words); cnt=Counter()
    for w in words:
        deg=sum(1 for i in range(n) if tuple(w[:i]+(1-w[i],)+w[i+1:]) in S)
        cnt[deg]+=1
    return cnt
x,y=sp.symbols('x y')
R=x**3*y+x**5*y**2/(1-x**2)
D=1-x*y-R+x*(y-1)*R**2
F=1+x*(1-y)+x**2*(1-y**2)-x*sp.diff(D,x)/D+x*(1-y)*(R+x*sp.diff(R,x))
NMAX=18
ser=sp.series(F,x,0,NMAX+1).removeO()
ok=True
for n in range(0,NMAX+1):
    c=sp.Poly(sp.expand(ser.coeff(x,n)),y) if n>0 else sp.Poly(sp.expand(ser.subs(x,0)),y)
    gf={m[0]:int(v) for m,v in c.terms()} if c.terms() else {}
    gf={k:v for k,v in gf.items() if v!=0}
    bf=dict(N(n)) if n>0 else {0:1}
    if gf!=bf: ok=False; print('MISMATCH',n,gf,bf)
    if n in (3,5,8,12,18): print(n,sorted(bf.items()))
print('all n<=%d match:'%NMAX, ok)
