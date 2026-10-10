import numpy as np, sys
from scipy.optimize import minimize, brentq
rng=np.random.default_rng(int(sys.argv[1]) if len(sys.argv)>1 else 0)
L=np.log(2)
def S(r):
    w=np.linalg.eigvalsh((r+r.conj().T)/2); w=w[w>1e-15]; return -np.sum(w*np.log(w))
def h(p):
    p=min(max(p,0.0),1.0); return 0.0 if p in (0.0,1.0) else -p*np.log(p)-(1-p)*np.log(1-p)
def hinv(y):
    if y<=0: return 0.0
    if y>=L: return 0.5
    return brentq(lambda p:h(p)-y,0,0.5)
def conv(a,b): return a*(1-b)+(1-a)*b
def bound(H1,H2):
    if H1+H2<=L: return h(conv(hinv(H1),hinv(H2)))
    return H1+H2-L+h(conv(hinv(L-H1),hinv(L-H2)))
def dm(v,d):
    A=(v[:d*d]+1j*v[d*d:2*d*d]).reshape(d,d); M=A@A.conj().T; return M/np.trace(M).real
def condH(rs,p):  # H(X|B) for cq state sum p_x |x><x| ⊗ rs[x]
    return sum(p[x]*0 for x in range(2)) + (sum(-p[x]*np.log(p[x]) for x in range(2) if p[x]>0) + sum(p[x]*S(rs[x]) for x in range(2))) - S(p[0]*rs[0]+p[1]*rs[1])
def gap(v,d):
    k=2*d*d
    r10,r11,r20,r21=dm(v[0:k],d),dm(v[k:2*k],d),dm(v[2*k:3*k],d),dm(v[3*k:4*k],d)
    p1=1/(1+np.exp(-v[4*k])); p2=1/(1+np.exp(-v[4*k+1]))
    P1=[1-p1,p1]; P2=[1-p2,p2]
    H1=condH([r10,r11],P1); H2=condH([r20,r21],P2)
    s0=P1[0]*P2[0]*np.kron(r10,r20)+P1[1]*P2[1]*np.kron(r11,r21)
    s1=P1[0]*P2[1]*np.kron(r10,r21)+P1[1]*P2[0]*np.kron(r11,r20)
    q0=np.trace(s0).real; q1=np.trace(s1).real
    HZ=(-q0*np.log(q0)-q1*np.log(q1)+q0*S(s0/q0)+q1*S(s1/q1))-S(s0+s1)
    return HZ-bound(max(H1,0),max(H2,0))
d=int(sys.argv[2]) if len(sys.argv)>2 else 2
best=1e9
for r in range(int(sys.argv[3]) if len(sys.argv)>3 else 40):
    v0=rng.standard_normal(8*d*d+2)*rng.choice([0.3,1,3])
    res=minimize(gap,v0,args=(d,),method='Nelder-Mead',options={'maxiter':20000,'xatol':1e-10,'fatol':1e-12})
    res=minimize(gap,res.x,args=(d,),method='BFGS')
    if res.fun<best: best=res.fun; print('restart',r,'min gap',repr(best),flush=True)
print('FINAL d=',d,'min gap',repr(best))
