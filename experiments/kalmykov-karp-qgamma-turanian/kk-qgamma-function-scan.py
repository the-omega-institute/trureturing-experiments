from mpmath import mp, mpf
mp.dps=25
def rgvals(s,q,M,N):
    qq=mpf(1)
    for j in range(N): qq*=(1-q**(j+1))
    p=mpf(1)
    for j in range(N): p*=(1-q**(s+j))
    r=p/(qq*(1-q)**(1-s)); out=[]; x=s
    for n in range(M):
        out.append(r)
        br=(1-q**x)/(1-q)
        if abs(x)<mpf('1e-20'):
            # pole: rg(x+1)=rg(1)... recompute directly at next point
            r=None
        else:
            r=r/br
        x+=1
        if r is None:
            p=mpf(1)
            for j in range(N): p*=(1-q**(x+j))
            r=p/(qq*(1-q)**(1-x))
    return out
def g(c,x):
    return sum(cn*x**n for n,cn in enumerate(c))
import sys
for q in [mpf(s) for s in sys.argv[1:]]:
    N=int(60/(1-q))+200; M=int(80/(1-q))+100; R=1/(1-q)
    for (mu,a,b) in [(-0.7,0.7,0.7),(-0.5,0.5,0.5),(-0.9,0.9,0.9),(-0.3,0.3,0.3),(-0.7,0.3,0.7),(-0.7,1.5,0.7)]:
        mu,a,b=mpf(mu),mpf(a),mpf(b)
        c0,ca,cb,cab=[rgvals(s,q,M,N) for s in (mu,mu+a,mu+b,mu+a+b)]
        worst=None
        for k in range(1,60):
            x=R*k/60
            A=g(ca,x)*g(cb,x); B=g(c0,x)*g(cab,x)
            r=(A-B)/(abs(A)+abs(B))
            if worst is None or r<worst[0]: worst=(r,x)
        print(float(q),float(mu),float(a),float(b),'min rel Delta',mp.nstr(worst[0],5),'x',mp.nstr(worst[1],4),'R',mp.nstr(R,4),flush=True)
