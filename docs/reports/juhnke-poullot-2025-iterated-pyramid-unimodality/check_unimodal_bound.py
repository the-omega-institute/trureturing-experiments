import random, itertools
from math import comb
def unimodal(a):
    a=list(a)
    while a and a[-1]==0: a.pop()
    while a and a[0]==0: a.pop(0)
    i=0
    while i+1<len(a) and a[i+1]>=a[i]: i+=1
    while i+1<len(a) and a[i+1]<=a[i]: i+=1
    return i==len(a)-1
def mulbin(q,m):
    out=[0]*(len(q)+m)
    for i,c in enumerate(q):
        if c:
            for j in range(m+1): out[i+j]+=c*comb(m,j)
    return out
bad=0; worst={}
random.seed(1)
for d in range(1,8):
    m=4*d*d
    trials=[]
    for _ in range(3000):
        q=[random.choice([0,0,1,random.randint(0,10**random.randint(0,12))]) for _ in range(d+1)]
        if q[d]==0: q[d]=random.randint(1,10**8)
        trials.append(q)
    for M in [10**k for k in range(0,30,3)]:
        for pos in range(d+1):
            q=[1]*(d+1); q[pos]=M; trials.append(q)
            q=[0]*(d+1); q[0]=1; q[pos]=M; q[d]=max(q[d],1); trials.append(q)
            q=[M]*(d+1); q[pos]=1; trials.append(q)
    for q in trials:
        if sum(q)==0: continue
        if not unimodal(mulbin(q,m)): bad+=1; print('FAIL d',d,'q',q)
    # minimal m that works for adversarial 1 + M X^d over all M
    mm=0
    for M in [10**k for k in range(0,40,2)]:
        q=[1]+[0]*(d-1)+[M]
        m2=0
        while not unimodal(mulbin(q,m2)): m2+=1
        mm=max(mm,m2)
    worst[d]=(mm,4*d*d)
print('bad',bad); print('min m for 1+M X^d (max over M) vs 4d^2:',worst)
