# Dekking–Keane arXiv:2202.13548v1 Conjecture 4: frequency of 1 in the fixed point x^(00) of kappa(ab) = a b (1-b) on pairs.
import math, sys
def word(N):
    x=[0,0]
    for i in range(2,N):
        q,r=divmod(i,3)
        x.append(x[2*q] if r==0 else (x[2*q+1] if r==1 else 1-x[2*q+1]))
    return x
N=int(sys.argv[1]) if len(sys.argv)>1 else 3**13
x=word(N)
# recurrence check x_{3n}=x_{2n}, x_{3n+1}=x_{2n+1}, x_{3n+2}=1-x_{2n+1}
for n in range(N//3):
    assert x[3*n]==x[2*n] and x[3*n+1]==x[2*n+1] and x[3*n+2]==1-x[2*n+1]
print('prefix', ''.join(map(str,x[:33])))
# coefficient functional d_k: sum(T^k v) = sum_j d_k(j) v_j, T(u,v)=(u,v,-v)
def T(v):
    out=[]
    for i in range(0,len(v),2): out+= [v[i],v[i+1],-v[i+1]]
    return out
for k in range(0,11):
    m=2**k
    d=[]
    for j in range(m):
        e=[0]*m; e[j]=1; w=e
        for _ in range(k): w=T(w)
        d.append(sum(w))
    E=sum(t*t for t in d); R=sum(d[r]*d[(r+1)%m] for r in range(m))
    exp_E = 1 if k==0 else 3**(k-1); exp_R = 1 if k==0 else 0
    assert E==exp_E and R==exp_R, (k,E,R)
    print(f'k={k} E_k={E} R_k={R}')
a=[1-2*t for t in x]; s=0; worst=0
for i,t in enumerate(a):
    s+=t
    if i+1>=1000: worst=max(worst,abs(s)/(i+1)**(math.log(math.sqrt(6))/math.log(3)))
print('freq of 1 at N=',N,':', sum(x)/N, ' max |S_N|/N^alpha (N>=1000) =', round(worst,4))
print('OK')
