# Asgarli-Yip arXiv:2110.07176 Conjecture 5.9: for odd prime p, a 2-dim F_p-subspace V of F_{p^4} with 1 in V
# is a clique in PP(p^4, p+1, I) for some I (|I|=(p+1)/2) iff V = F_p + aF_p with a = g^{(p+1)k}, k odd.
# Clique for some I  <=>  V\{0} meets at most (p+1)/2 of the p+1 cyclotomic classes C_j = g^j <g^{p+1}>.
import itertools, sys
def run(p):
    # GF(p^4) = F_p[x]/(f), elements as tuples of 4 coefficients
    def mul(a,b,f):
        r=[0]*7
        for i in range(4):
            for j in range(4): r[i+j]=(r[i+j]+a[i]*b[j])%p
        for k in range(6,3,-1):  # reduce x^k, f monic: x^4 = -(f0+f1x+f2x^2+f3x^3)
            c=r[k]
            if c:
                r[k]=0
                for i in range(4): r[k-4+i]=(r[k-4+i]-c*f[i])%p
        return tuple(r[:4])
    one=(1,0,0,0); Q=p**4-1
    def power(a,e,f):
        res=one; base=a
        while e:
            if e&1: res=mul(res,base,f)
            base=mul(base,base,f); e>>=1
        return res
    primes=[q for q in range(2,Q+1) if Q%q==0 and all(q%d for d in range(2,int(q**.5)+1))]
    for f in itertools.product(range(p),repeat=4):
        x=(0,1,0,0)
        if f[0]==0: continue
        if power(x,Q,f)!=one: continue
        if all(power(x,Q//q,f)!=one for q in primes): break
    g=(0,1,0,0)
    # discrete log table
    log={}; e=one
    for i in range(Q):
        log[e]=i; e=mul(e,g,f)
    assert len(log)==Q
    def add(a,b): return tuple((x+y)%p for x,y in zip(a,b))
    def smul(c,a): return tuple((c*x)%p for x in a)
    cls=lambda v: log[v]%(p+1)
    # predicted set: subspaces F_p + aF_p with a = g^{(p+1)k}, k odd; identify a subspace by its sorted element set
    def span(b): return frozenset(add(smul(s,one),smul(t,b)) for s in range(p) for t in range(p))
    elems=[v for v in log]
    Fp=set(smul(c,one) for c in range(p))
    subspaces={}
    for b in elems:
        if b in Fp: continue
        V=span(b)
        if V not in subspaces: subspaces[V]=b
    pred=set()
    gp1=power(g,p+1,f)
    for k in range(1,Q//(p+1)+1,2):
        a=power(gp1,k,f)
        if a in Fp: continue
        pred.add(span(a))
    actual=set(V for V in subspaces if len(set(cls(v) for v in V if any(v)))<=(p+1)//2)
    return len(subspaces), len(actual), len(pred), actual==pred
for p in [int(a) for a in sys.argv[1:]] or [3,5,7]:
    print("p",p,"subspaces_with_1",*run(p),flush=True)
