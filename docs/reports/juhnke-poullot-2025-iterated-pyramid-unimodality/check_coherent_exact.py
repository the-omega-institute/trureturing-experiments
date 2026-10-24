"""Coherent/monotone path counts on iterated pyramids via HiGHS LPs (scipy).
edge [u,v]: exists c with <c,u>=<c,v> >= <c,w>+1 for all other vertices w.
coherent path v0..vl: exists w with dc_i <w,u-v_i> - <c,u-v_i><w,v_{i+1}-v_i> <= -1 for all steps i, u not in {v_i,v_{i+1}}."""
import numpy as np, itertools, random, sys
from collections import Counter
from scipy.optimize import linprog
def feasible(A_ub,b_ub,A_eq=None,b_eq=None,n=None):
    r=linprog(np.zeros(n),A_ub=A_ub,b_ub=b_ub,A_eq=A_eq,b_eq=b_eq,bounds=[(-1e6,1e6)]*n,method='highs')
    return r.status==0
def edges(V):
    V=np.array(V,float); n,d=V.shape; E=[]
    for i,j in itertools.combinations(range(n),2):
        others=[k for k in range(n) if k not in (i,j)]
        A=np.array([V[k]-V[i] for k in others]); b=-np.ones(len(others))
        if feasible(A,b,np.array([V[i]-V[j]]),np.zeros(1),d): E.append((i,j))
    return E
def paths(V,c,E):
    h=np.array(V,float)@np.array(c,float); adj={i:[] for i in range(len(V))}
    for i,j in E:
        assert abs(h[i]-h[j])>1e-9
        a,b=(i,j) if h[i]<h[j] else (j,i); adj[a].append(b)
    s=int(np.argmin(h)); t=int(np.argmax(h)); out=[]
    def go(p):
        if p[-1]==t: out.append(tuple(p)); return
        for y in adj[p[-1]]: go(p+[y])
    go([s]); return out
def coherent(V,c,p):
    V=np.array(V,float); c=np.array(c,float); n,d=V.shape; rows=[]
    for a,b in zip(p,p[1:]):
        dc=c@(V[b]-V[a])
        for u in range(n):
            if u in (a,b): continue
            du=V[u]-V[a]; rows.append(dc*du-(c@du)*(V[b]-V[a]))
    A=np.array(rows); return feasible(A,-np.ones(len(rows)),n=d)
def counts(V,c):
    E=edges(V); P=paths(V,c,E); L=max(len(q)-1 for q in P)
    m=Counter(len(q)-1 for q in P); h=Counter(len(q)-1 for q in P if coherent(V,c,q))
    return [m.get(l,0) for l in range(L+1)],[h.get(l,0) for l in range(L+1)]
def pyr(V,c):
    mn=min(np.dot(v,c) for v in V)
    return [tuple([0.0]*(len(c)+1))]+[tuple(v)+(1.0,) for v in V], tuple(c)+(1.0-mn,)
def vertices(V):
    E=edges(V); keep=sorted(set(i for e in E for i in e)); return [V[i] for i in keep]
def mul(a,m):
    out=list(a)
    for _ in range(m): out=[x+y for x,y in zip(out+[0],[0]+out)]
    return out
def unimodal(a):
    i=0
    while i+1<len(a) and a[i+1]>=a[i]: i+=1
    while i+1<len(a) and a[i+1]<=a[i]: i+=1
    return i==len(a)-1
random.seed(3)
rand4=[tuple(random.randint(-20,20)/random.randint(1,5) for _ in range(4)) for _ in range(9)]
random.seed(11)
rand4b=[tuple(random.gauss(0,1) for _ in range(4)) for _ in range(10)]
rand5=[tuple(random.gauss(0,1) for _ in range(5)) for _ in range(10)]
ok=True
for name,V0,c0,K in [('rand4',rand4,(1,2,5,11),3),('rand4b',rand4b,(1,2.3,5.7,11.1),3),('rand5',rand5,(1,2.3,5.7,11.1,23.9),2)]:
    V=vertices(V0); c=tuple(float(x) for x in c0); base=None
    for k in range(K+1):
        if k: V,c=pyr(V,c)
        mono,coh=counts(V,c)
        print(f'{name} k={k}: |V|={len(V)} mono={mono} coh={coh} strict={mono!=coh} unimodal=({unimodal(mono)},{unimodal(coh)})',flush=True)
        if k==1: base=(mono,coh)
        if k>=2:
            pm,pc=mul(base[0],k-1),mul(base[1],k-1)
            print(f'  product X(1+X)^(k-1)Q: mono {mono==pm} coh {coh==pc}',flush=True); ok&=mono==pm and coh==pc
print('HI_OK' if ok else 'HI_FAIL')
