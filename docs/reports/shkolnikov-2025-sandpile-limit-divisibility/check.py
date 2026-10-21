"""Check for Shkolnikov's expectation (arXiv:2506.09912v3, Section 3) that the direct limit of rectangular sandpile
groups is divisible.

G(p, q) is represented as in Remark 1 by Q/Z-valued functions on Gamma(p, q) = (0,p) x (0,q) cap Z^2, extended by 0,
with 4 f(x,y) = sum of the four neighbours; elements are produced as Delta^-1 v mod 1 for integer v.  For a given
f and k, the script builds g on Gamma(p, qP) by the row recurrence r_{j+1} = B r_j - r_{j-1} (B = 4I - A), starting
from s_0 = 0 and s_1 = r_1/k, with P the period of the transfer map, and checks exactly that g is harmonic and that
k g equals the source's signed reflected-tiling embedding of f.  It also reproduces orders from the source's Table 1.
"""
from fractions import Fraction as F
import itertools, random
random.seed(7)
def solve(A,b):
    n=len(A); M=[row[:]+[b[i]] for i,row in enumerate(A)]
    for c in range(n):
        piv=next(r for r in range(c,n) if M[r][c]!=0); M[c],M[piv]=M[piv],M[c]
        iv=1/M[c][c]; M[c]=[v*iv for v in M[c]]
        for r in range(n):
            if r!=c and M[r][c]!=0:
                f=M[r][c]; M[r]=[M[r][j]-f*M[c][j] for j in range(n+1)]
    return [M[i][n] for i in range(n)]
def sites(p,q): return [(x,y) for x in range(1,p) for y in range(1,q)]
def lap(p,q):
    S=sites(p,q); idx={s:i for i,s in enumerate(S)}
    A=[[F(0)]*len(S) for _ in S]
    for (x,y),i in idx.items():
        A[i][i]=F(4)
        for d in ((1,0),(-1,0),(0,1),(0,-1)):
            t=(x+d[0],y+d[1])
            if t in idx: A[i][idx[t]]=F(-1)
    return S,idx,A
def frac(v): return v-(v.numerator//v.denominator)
def harmonic(p,q,f):
    for (x,y) in sites(p,q):
        s=4*f.get((x,y),F(0))-sum(f.get((x+a,y+b),F(0)) for a,b in ((1,0),(-1,0),(0,1),(0,-1)))
        if frac(s)!=0: return False
    return all(0<x<p and 0<y<q for (x,y) in f)
def random_element(p,q):
    S,idx,A=lap(p,q); v=[F(random.randint(-3,3)) for _ in S]
    sol=solve(A,v); return {s:frac(sol[idx[s]]) for s in S}
def unfold(p,x):
    k=x//p; return x-k*p if k%2==0 else (k+1)*p-x
def embed(p,q,m,n,f):
    g={}
    for x in range(1,m*p):
        for y in range(1,n*q):
            v=f.get((unfold(p,x),unfold(q,y)),F(0)); s=(-1)**(x//p+y//q)
            g[(x,y)]=frac(s*v)
    return g
def root(p,q,f,k):
    # rows r_j over j=0..q (r_0=r_q=0), as vectors over x=1..p-1; choose s_1 with k s_1 = r_1, iterate until return
    r1=[f.get((x,1),F(0)) for x in range(1,p)]
    s_prev=[F(0)]*(p-1); s_cur=[frac(v/k) for v in r1]
    rows=[s_prev,s_cur]; start=(tuple(s_prev),tuple(s_cur))
    def step(u,v):
        w=[frac(4*v[i]-(v[i-1] if i>0 else 0)-(v[i+1] if i<p-2 else 0)-u[i]) for i in range(p-1)]
        return w
    P=None
    for j in range(1,100000):
        w=step(rows[-2],rows[-1]); rows.append(w)
        if (tuple(rows[-2]),tuple(rows[-1]))==start: P=len(rows)-2; break
    # P = period of transfer from state (s_0,s_1); need qP rows: s_{qP}=0
    n=P  # use n=P tiles vertically: height q*P
    H=q*n; s=rows
    while len(s)<H+2: s.append(step(s[-2],s[-1]))
    assert all(v==0 for v in s[H]), "s_qP nonzero"
    g={(x,y):s[y][x-1] for x in range(1,p) for y in range(1,H)}
    return n,g
ok=True; cnt=0
for (p,q) in [(2,2),(3,2),(2,3),(3,3),(4,3)]:
    for k in (2,3,5):
        for _ in range(3):
            f=random_element(p,q); assert harmonic(p,q,f)
            n,g=root(p,q,f,k)
            e=embed(p,q,1,n,f)
            okh=harmonic(p,q*n,g) and harmonic(p,q*n,e)
            eq=all(frac(k*g.get(z,F(0))-e.get(z,F(0)))==0 for z in set(g)|set(e))
            ok&=okh and eq; cnt+=1
print("cases", cnt, "all roots harmonic and k*g = embed(f):", ok)
def order(p, q):
    S, idx, A = lap(p, q)
    n = len(S); M = [row[:] for row in A]; det = F(1)
    for c in range(n):
        piv = next(r for r in range(c, n) if M[r][c] != 0)
        if piv != c: M[c], M[piv] = M[piv], M[c]; det = -det
        det *= M[c][c]
        for r in range(c + 1, n):
            f_ = M[r][c] / M[c][c]
            M[r] = [M[r][j] - f_ * M[c][j] for j in range(n)]
    return det
tab = {(2, 2): 4, (2, 3): 15, (3, 3): 192, (2, 4): 56, (3, 4): 2415, (4, 4): 100352}
for pq, v in tab.items():
    o = order(*pq)
    print(f"|G{pq}| = {o} (source Table 1: {v})")
    ok &= o == v
print("ALL_CHECKS", ok)
raise SystemExit(0 if ok else 1)
