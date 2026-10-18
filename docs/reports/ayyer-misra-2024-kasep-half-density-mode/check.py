"""Ayyer-Misra arXiv:2310.03343v2 Section 4.5: exact numerator polynomials of pi(eta_1=k) and pi(eta_1=k+1) at
n=(2k+1)L/2, L<=10, K<=6. Checks the first nonzero coefficient of their difference: offset min(k+1,K-k) above the
lowest degree of Z, value sign(K-2k-1)*3*C(L-1,L/2-2), Z leading coefficient C(L,L/2); tie when 2k+1=K.
Also prints the source's claim-1 counterexample above half filling (L=4, K=2, k=1, n=6).
"""
import sys
import sympy as sp
t=sp.symbols('t')
def qbin(K,a):
    if not 0<=a<=K: return sp.Integer(0)
    num=sp.prod([1-t**(K-i) for i in range(a)]); den=sp.prod([1-t**(i+1) for i in range(a)])
    q,r=sp.div(sp.Poly(sp.expand(num),t),sp.Poly(sp.expand(den),t)); assert r.is_zero
    return q.as_expr()
def weights(K): return {e: sp.expand(t**(((e-1)*(e-2))//2)*qbin(K,e)) for e in range(K+1)}
def tail(L,K,s):
    w=weights(K); dp={0:sp.Integer(1)}
    for _ in range(L): 
        nd={}
        for a,v in dp.items():
            for e in range(K+1):
                if a+e<=s: nd[a+e]=sp.expand(nd.get(a+e,0)+v*w[e])
        dp=nd
    return dp.get(s,sp.Integer(0))
def lowest(p):
    P=sp.Poly(sp.expand(p),t)
    if P.is_zero: return None,0
    m=min(P.monoms())[0]; return m, P.coeff_monomial(t**m)
ok=True; rows=0
for L in (4,6,8,10):
    r=L//2
    for K in range(1,7):
        w=weights(K)
        for k in range(K):
            n=k*L+r
            rest_k=tail(L-1,K,n-k); rest_k1=tail(L-1,K,n-k-1)
            A=sp.expand(w[k]*rest_k); B=sp.expand(w[k+1]*rest_k1); Z=tail(L,K,n)
            eZ,cZ=lowest(Z); eD,cD=lowest(A-B)
            rows+=1
            if 2*k+1==K:
                good = eD is None
            else:
                s=min(k+1,K-k); sign=1 if K>2*k+1 else -1
                good = (eD is not None and eD-eZ==s and cD==sign*3*sp.binomial(L-1,r-2) and cZ==sp.binomial(L,r))
            if not good:
                ok=False; print('MISMATCH',L,K,k,'eD-eZ=',None if eD is None else eD-eZ,'cD=',cD,'cZ=',cZ)
print('rows',rows,'GPT_PRO_FORMULA_OK',ok)

w=weights(2); A=sp.expand(w[1]*tail(3,2,5)); B=sp.expand(w[2]*tail(3,2,4))
print("L=4 K=2 k=1 n=6 numerators:", sp.Poly(A,t).all_coeffs()[::-1], sp.Poly(B,t).all_coeffs()[::-1])
sys.exit(0 if ok else 1)
