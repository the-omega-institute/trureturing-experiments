import itertools, sympy as sp, sys
def valid(w):
    n=len(w)
    if n==1: return w==(0,)
    if 1 not in w: return True
    if 0 not in w: return False
    i=next(i for i in range(n) if w[i]==1 and w[i-1]==0)
    r=w[i:]+w[:i]; j=0
    while j<n:
        a=0
        while j<n and r[j]==1: a+=1; j+=1
        b=0
        while j<n and r[j]==0: b+=1; j+=1
        if b<=a: return False
    return True
def cube_poly_subcube(n):
    V={w for w in itertools.product((0,1),repeat=n) if valid(w)}
    c={}
    for v in V:
        ones=[i for i in range(n) if v[i]]
        for k in range(len(ones)+1):
            for S in itertools.combinations(ones,k):
                ok=True
                for m in range(1,2**k):
                    u=list(v)
                    for t,s in enumerate(S):
                        if m>>t&1: u[s]=0
                    if tuple(u) not in V: ok=False;break
                if ok: c[k]=c.get(k,0)+1
    return V,c
def graph_counts(V):
    V=list(V); idx={v:i for i,v in enumerate(V)}; n=len(V[0])
    adj=[set() for _ in V]
    for v in V:
        for i in range(n):
            u=list(v); u[i]^=1; u=tuple(u)
            if u in idx: adj[idx[v]].add(idx[u])
    E=sum(len(a) for a in adj)//2
    # 4-cycles: count pairs of vertices at distance 2 with >=2 common neighbours, each square counted twice
    sq=0
    for a in range(len(V)):
        for b in range(a+1,len(V)):
            cm=len(adj[a]&adj[b])
            sq+=cm*(cm-1)//2
    return E, sq//2
x,z=sp.symbols('x z')
D=1-z-z**2-x*z**3-x*(1+x)*z**5
G=(z+2*z**2+3*x*z**3+5*x*(1+x)*z**5)/D-2*z**2/(1-z**2)
N=int(sys.argv[1]) if len(sys.argv)>1 else 12
ser=sp.series(G,z,0,N+1).removeO()
ok=True
for n in range(1,N+1):
    V,c=cube_poly_subcube(n)
    poly=sp.expand(sum(cnt*x**k for k,cnt in c.items()))
    gf=sp.expand(ser.coeff(z,n))
    E,sq=graph_counts(V)
    good=(sp.expand(poly-gf)==0) and c.get(0,0)==len(V) and c.get(1,0)==E and c.get(2,0)==sq
    ok&=good
    print(n,len(V),poly,'| GF',gf,'| edges',E,'squares',sq,'OK' if good else 'MISMATCH')
print('ALL_OK' if ok else 'FAIL'); sys.exit(0 if ok else 1)
