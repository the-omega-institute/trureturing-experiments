import itertools, random
from fractions import Fraction
random.seed(11)
def span(gens,n):
    S={tuple([0]*n)}
    for g in gens: S|={tuple((a+b)%2 for a,b in zip(s,g)) for s in S}
    return S
def rank2(vs,n):
    rows=[list(v) for v in vs]; r=0
    for c in range(n):
        p=next((i for i in range(r,len(rows)) if rows[i][c]),None)
        if p is None: continue
        rows[r],rows[p]=rows[p],rows[r]
        for i in range(len(rows)):
            if i!=r and rows[i][c]: rows[i]=[(a+b)%2 for a,b in zip(rows[i],rows[r])]
        r+=1
    return r
def dual(S,n): return [v for v in itertools.product([0,1],repeat=n) if all(sum(a*b for a,b in zip(v,s))%2==0 for s in S)]
def poly_from_exps(exps,n):
    c=[0]*(n+1)
    for e in exps: c[e]+=1
    return c
def mult(c,a):  # root multiplicity of integer polynomial c (coeff list low->high) at integer a, by synthetic division
    c=[Fraction(x) for x in c]; m=0
    while any(c):
        # evaluate
        r=Fraction(0); q=[]
        for coef in reversed(c): r=r*a+coef; q.append(r)
        if r!=0: break
        q=q[:-1]; q.reverse(); c=q; m+=1
    return m
def analyse(GX,GZ,n):
    one=tuple([1]*n); GXs=set(GX); GZs=set(GZ)
    GXp=dual(GZ,n); GZp=dual(GX,n)
    d=min(sum(1 for i in range(n) if x[i] or zz[i]) for x in GXp for zz in GZp if not (x in GXs and zz in GZs))
    P=poly_from_exps([n-sum(g) for g in GX],n); Q=poly_from_exps([sum(g) for g in GX],n)
    ev=lambda c,a: sum(Fraction(x)*Fraction(a)**i for i,x in enumerate(c))
    res={}
    for a in [0,1,-1]:
        assert ev(Q,a)!=0 and ev(P,a)==a*ev(Q,a), ('fix',a)
        res[a]=mult([p-a*q for p,q in zip(P,Q)],a)
    Qr=list(reversed(Q)); Pr=list(reversed(P))   # degree-n reversals
    assert Pr[0]!=0
    res['inf']=mult(Qr,0)
    dX=min(n-sum(g) for g in GX); dZ=min(n-sum(h) for h in GZ)
    return d,res,dX,dZ
def steane():
    H=[(1,0,1,0,1,0,1),(0,1,1,0,0,1,1),(0,0,0,1,1,1,1)]
    S=list(span(H,7)); return S,S,7
def random_css(n):
    for _ in range(5000):
        kx=random.randint(1,(n-1)//2)
        gens=[tuple(random.randint(0,1) for _ in range(n)) for _ in range(kx)]
        gens=[g if sum(g)%2==0 else tuple(list(g[:-1])+[1-g[-1]]) for g in gens]
        if rank2(gens,n)<kx: continue
        GX=span(gens,n); one=tuple([1]*n)
        if one in GX: continue
        cand=[v for v in dual(list(GX),n) if sum(v)%2==0 and v!=one]
        for _ in range(100):
            if len(cand)<n-1-kx: break
            gz=random.sample(cand,n-1-kx)
            if rank2(gz,n)==n-1-kx:
                GZ=span(gz,n)
                if one not in GZ: return list(GX),list(GZ),n
    return None
out=[]
GX,GZ,n=steane(); out.append(('steane',n)+analyse(GX,GZ,n))
cnt=0
while cnt<80:
    n=random.choice([5,7,9,11])
    r=random_css(n)
    if r is None: continue
    GX,GZ,n=r; out.append(('rand',n)+analyse(GX,GZ,n)); cnt+=1
bad=[o for o in out if min(o[3].values())<o[2]]
exact=[o for o in out if o[3][0]==o[4] and o[3]['inf']==o[4] and o[3][1]==o[5] and o[3][-1]==o[5]]
print('codes',len(out),'violations',len(bad),'exact (dX at 0,inf; dZ at +-1)',len(exact))
print('steane',out[0])
from collections import Counter
print('d values',Counter(o[2] for o in out))
print('min over codes of (order - d)', min(min(o[3].values())-o[2] for o in out))
