"""Exact check of the decomposition H_n = (C_2n x C_2^(n-2)) x| C_n for the Cui-Mehta-Mousavi-Nezhadi groups.
Elements of H_n are monomial matrices diag(q^e_0..q^e_{n-1}) X^j with q = exp(i*pi/n) (order 2n); stored as (e mod 2n tuple, j mod n).
X e_i = e_{i+1}; A0 = X, A1 = q D_0 X where D_0 = diag(-1 at 0) = q^n at coordinate 0."""
import sys
from itertools import product
def mul(a,b,n):
    (ea,ja),(eb,jb)=a,b
    # (Da X^ja)(Db X^jb) = Da (X^ja Db X^-ja) X^(ja+jb); X^j diag(d_k) X^-j = diag(d_{k-j})
    e=tuple((ea[k]+eb[(k-ja)%n])%(2*n) for k in range(n))
    return (e,(ja+jb)%n)
def ident(n): return (tuple([0]*n),0)
def gen_group(gens,n):
    G={ident(n)}; frontier=[ident(n)]
    while frontier:
        nf=[]
        for g in frontier:
            for h in gens:
                x=mul(g,h,n)
                if x not in G: G.add(x); nf.append(x)
        frontier=nf
    return G
def power(a,k,n):
    r=ident(n)
    for _ in range(k): r=mul(r,a,n)
    return r
def inv(a,n):
    e,j=a
    # inverse of D X^j is X^-j D^-1 = (X^-j D^-1 X^j) X^-j
    e2=tuple((-e[(k+j)%n])%(2*n) for k in range(n))
    return (e2,(-j)%n)
ok=True
for n in range(2,10):
    X=(tuple([0]*n),1)
    A1=(tuple([1+n if k==0 else 1 for k in range(n)]),1)   # q*D_0 = q^(1+n) at coord 0, q elsewhere
    H=gen_group([X,A1],n)
    order_ok = len(H)==n*n*2**(n-1)
    s=mul(A1,inv(X,n),n)                       # s_0 = A1 X^-1 = q D_0
    sI=[mul(mul(power(X,i,n),s,n),inv(power(X,i,n),n),n) for i in range(n)]
    t=[mul(sI[i],inv(s,n),n) for i in range(n)]   # t_i = s_i s^-1, t_0 = 1
    # t_{n-1} = s^n prod_{i=1}^{n-2} t_i
    prod_t=ident(n)
    for i in range(1,n-1): prod_t=mul(prod_t,t[i],n)
    tlast_ok = t[n-1]==mul(power(s,n,n),prod_t,n)
    # Phi: C_2n x C_2^(n-2) -> K injective and onto K = <s_i>
    K=gen_group(sI,n)
    imgs=set()
    for a in range(2*n):
        for b in product(range(2),repeat=n-2):
            g=power(s,a,n)
            for i,bi in enumerate(b):
                if bi: g=mul(g,t[i+1],n)
            imgs.add(g)
    phi_ok = len(imgs)==2*n*2**(n-2) and imgs==K
    # action: X s X^-1 = s t_1 (t_1 := t[1] if n>=3 else s^2), X t_i X^-1 = t_{i+1} t_1
    t1 = t[1] if n>=3 else power(s,2,n)
    conj=lambda g: mul(mul(X,g,n),inv(X,n),n)
    act_ok = conj(s)==mul(s,t1,n) and all(conj(t[i])==mul(t[i+1],t1,n) for i in range(1,n-1))
    # trivial intersection K cap <X>
    cyc={power(X,j,n) for j in range(n)}
    inter_ok = (K & cyc)=={ident(n)}
    print(n,len(H),order_ok,tlast_ok,phi_ok,act_ok,inter_ok)
    ok = ok and order_ok and tlast_ok and phi_ok and act_ok and inter_ok
print('ALL_OK' if ok else 'FAIL'); sys.exit(0 if ok else 1)
