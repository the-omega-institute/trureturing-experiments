# Outward-rounded rational certificate (standard library only).
# Kalmykov-Karp, arXiv:1702.03788, Theorem 3 / Conjecture 1, function clause:
# Delta_g(alpha,beta;x) = g(mu+alpha;x) g(mu+beta;x) - g(mu;x) g(mu+alpha+beta;x) > 0 for 0 < x < R_g.
# Witness: g_n = 1, q = r^10 with r = 93/100, mu = -7/10, alpha = 3/10, beta = 7/10, x = X.
# g(s;x) = sum_n x^n rg(s+n), rg(s) = (q^s;q)_inf / ((q;q)_inf (1-q)^(1-s)), rg(s+1) = rg(s)/[s]_q.
# With s in {-7,-4,0,3}/10 the (1-q)^(1-s) factors are irrational, so work with
# h(s;x) = (q;q)_inf * (1-q)^(1-s) * g(s;x) = sum_n (q^s;q)_inf * x^n / prod_{j<n} [s+j]_q,
# and note (1-q)^{(1-s1)+(1-s2)} is the same for both products: s(mu+a)+s(mu+b) = s(mu)+s(mu+a+b).
# So sign(Delta) = sign(h(mu+a) h(mu+b) - h(mu) h(mu+a+b)).  For s = 0 (pole), rg(0) = 0 and the series starts at n = 1
# with rg(1) = 1/(q;q)_inf... handled exactly: (q^0;q)_inf = 0, but (q^0;q)_inf/[0]_q -> (q;q)_inf (1-q) is the limit;
# so for s = 0 use h(0;x) = sum_{n>=1} x^n (1-q) (q;q)_inf / prod_{1<=j<n} [j]_q ... and (1-q)^(1-0) factor as above.
from fractions import Fraction as F
import sys
r=F(93,100); q=r**10; D=1<<160
def dn(v): return F(v.numerator*D//v.denominator, D)
def up(v): return F(-((-v.numerator*D)//v.denominator), D)
def qpow(t): return r**t            # q^(t/10)
def br(t): return (1-qpow(t))/(1-q) # [t/10]_q
def poch_bounds(t,N=120):           # bounds for (q^(t/10);q)_inf, t>0
    lo=F(1); hi=F(1); a=qpow(t)
    xl=a; xh=a
    for _ in range(N):
        lo=dn(lo*(1-xh)); hi=up(hi*(1-xl)); xl=dn(xl*q); xh=up(xh*q)
    return dn(lo*(1-xh/(1-q))), hi
def series_bounds(t,X,P,K=160):
    # sum_{n>=0} X^n / prod_{j<n}[ (t+10j)/10 ]_q times P ; P = (lo,hi) of the prefactor
    # terms: c_0 = 1, c_{n+1} = c_n * X / [(t+10n)/10]_q ; tail after K with ratio bound X/[(t+10K)/10]_q <= X/[(t+10K)/10]
    s=F(0); c=F(1); out_lo=None
    for n in range(K):
        s+=c; c=c*X/br(t+10*n)
        c=up(c) if c>0 else -up(-c)   # magnitude rounded up; error per term <= 1/D
    # from n=K on, all terms share the sign of c and decay at ratio <= rho < 1 ([.]_q is increasing)
    rho=X/br(t+10*K)
    assert 0<rho<1 and t+10*K>0
    tail=c/(1-rho)                  # the true tail lies between 0 and this value
    fudge=F(K,D)                    # accumulated outward rounding of the computed terms
    S_lo=min(s,s+tail)-fudge; S_hi=max(s,s+tail)+fudge
    lo,hi=P
    return (min(lo*S_lo,hi*S_lo), max(lo*S_hi,hi*S_hi))
X=F(int(sys.argv[1]) if len(sys.argv)>1 else 93,100)
# s = mu+a = -4/10, mu+b = 0, mu = -7/10, mu+a+b = 3/10.
# (q^s;q)_inf for s<0: (q^s;q)_inf = (1-q^s) * (q^(s+1);q)_inf, with 1-q^s < 0.
def pref(t):   # bounds for (q^(t/10);q)_inf, any t not in -10N
    if t>0: return poch_bounds(t)
    lo,hi=poch_bounds(t+10); f=1-qpow(t)
    return (min(f*lo,f*hi),max(f*lo,f*hi))
hA=series_bounds(-4,X,pref(-4))
# s = 0: rg(0)=0; h(0;x) = (q;q)_inf (1-q) sum_{n>=1} x^n / prod_{1<=j<n} [j]_q  (limit of (q^s;q)_inf/[s]_q)
qq=poch_bounds(10)
tmp=series_bounds(10,X,qq)          # sum_{m>=0} X^m/prod_{j<m}[1+j]  * (q;q)_inf
hB=(tmp[0]*X*(1-q), tmp[1]*X*(1-q))
hC=series_bounds(-7,X,pref(-7))
hD=series_bounds(3,X,pref(3))
def mul(a,b):
    v=[a[0]*b[0],a[0]*b[1],a[1]*b[0],a[1]*b[1]]; return (min(v),max(v))
A=mul(hA,hB); B=mul(hC,hD)
print('X',X,'A in',float(A[0]),float(A[1]),' B in',float(B[0]),float(B[1]))
print('certified Delta<0:', A[1] < B[0])
