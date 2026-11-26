
from fractions import Fraction
from flint import arb, ctx

ctx.prec = 160
def q(value):
    fraction = Fraction(value)
    return arb(fraction.numerator)/fraction.denominator

def matrix(rows):
    return [[q(i) for i in row] for row in rows]

beta = list(map(q, ['41/10','-19/6',-7]))
B = matrix([['199/50','27/10','44/5'],['9/10','35/6',12],['11/10','9/2',-26]])
D = matrix([['17/10','1/2','3/2'],['3/2','3/2','1/2'],[2,2,0]])
C = matrix([['17/20','9/4',8],['1/4','9/4',8],['9/4','9/4',0]])
# (4,2,1) decomposes as a colored doublet plus a colorless doublet.
indices = matrix([['1/10','3/2',1],['3/10','1/2',0]])
casimir = matrix([['1/60','3/4','4/3'],['3/20','3/4',0]])
d = [sum(row[i] for row in indices)/3 for i in range(3)]
BD = [[sum(row[i]*(4*c[j]+q('2/3')*[0,2,3][i]*int(i==j))
           for row,c in zip(indices,casimir)) for j in range(3)] for i in range(3)]


from math import ceil
loop=16*arb.pi()**2
G0=[1/q(x)**2 for x in ('.467453','.63811','1.0583')]
yc=['.00000615','.00311','.8616','.0000135','.000268','.01401','.0000028227','.00059465','.010101']
y0=[2*q(x).log() for x in yc]
for k,err in [(3,'.0000002'),(4,'.000003')]:
    y0[k]=2*arb(yc[k],err).log()
def positive(x):
    if x<0:return arb(0)
    if x>0:return x
    return arb(0).union(x.upper())
def exact_exposure(t,a,offsets):
    u=t-a
    if u>q('.55'):return 9*u,arb(9)
    if u<q('-.768'):return arb(0),arb(0)
    if offsets is not None:
        return sum(positive(u-s) for s in offsets),sum(1 if u>s else 0 if u<s else arb(0).union(1) for s in offsets)
    A0=q('-.768'); B0=q('.55')
    theta=B0/(B0-A0)
    ub=9*(theta*positive(u-A0)+(1-theta)*positive(u-B0))
    return (9*positive(u)).lower().union(ub.upper()),arb(0).union(9)
gref=[q('.29'),q('.30'),q('.36')]
eff=[s+sum(z*g for z,g in zip(row,gref))/loop for s,row in zip(d,BD)]
vflow=[b+9*s for b,s in zip(beta,eff)]
rho=(vflow[1]-vflow[2])/(vflow[0]-vflow[2])
scalar_r=eff[1]-rho*eff[0]-(1-rho)*eff[2]
BF=[[6*b for b in row] for row in matrix([['19/15','3/5','44/15'],['1/5','49/3','4'],['11/30','3/2','76/3']])]

from fractions import Fraction as F

from functools import cache
gf=[q('.33')]*3
hf=[8+sum(b*g for b,g in zip(row,gf))/loop for row in BF]

def rates(gs,y,n):
    ys=[v.exp() for v in y]
    tr=[sum(ys[i:i+3]) for i in (0,3,6)]
    T=3*(tr[0]+tr[1])+tr[2]
    yp=[]
    for i,row in enumerate(C):
        opp=arb(0).union(tr[1-i]) if i<2 else arb(0)
        gy=sum(v*g for v,g in zip(row,gs))
        yp.extend((3*(v-opp)+2*(T-gy))/loop for v in ys[i*3:i*3+3])
    xp=[-2*(sum(B[i][j]*gs[j]+n*BD[i][j]*(gs[j]-gref[j])-D[i][j]*tr[j] for j in range(3)))/loop**2 for i in range(3)]
    fp=[-2*sum(v*(g-r) for v,g,r in zip(row,gs,gf))/loop**2 for row in BF]
    return xp,yp,fp

def physical(X,t,E,Z):
    return [v-2*(b*t+e*E+f*Z)/loop for v,b,e,f in zip(X,beta,eff,hf)]
def couplings(X,t,E,Z):
    G=physical(X,t,E,Z)
    assert all(v>1/(4*arb.pi()) for v in G)
    return [1/v for v in G]


raw_logs="""2776802599 2777383826
3291478277 3292006854
4011892601 4012355958
8044699025 8044928062
8045680522 8045909522
8046082531 8046311516
10458774531 10458929781
10536799409 10536952791
11882531321 11882656410
61116210131 61116210787
61116210131 61116210787
61116210131 61116210788
61116308348 61116309005
61116308390 61116309046
61116308396 61116309052
61401258324 61401258962"""

from fractions import Fraction
logs=[q(Fraction(int(a),10**10)).union(q(Fraction(int(b),10**10))) for a,b in
      (line.split() for line in raw_logs.strip().splitlines())]
def exposure(t,s,ls):
    v=t-s
    E=arb(0);n=arb(0)
    for w in ls:
        if v<w:break
        E+=positive(v-w)
        n+=1 if v>w else arb(0).union(1)
    return E,n
def run(al,au,tl,tu,steps=4500):
    a=q(al).union(q(au))
    tau=q(tl).union(q(tu))
    # Exactly 15 active at the crossing. u=tau-s=log(Q/M).
    W=9*(tau-a)
    u=(W+sum(logs[:15]))/15
    assert logs[14]<u and u<logs[15]
    X=G0[:]
    y=y0[:]
    dt=tau/steps
    for i in range(steps):
        fraction=q(Fraction(i-steps,steps)).union(q(Fraction(i+1-steps,steps)))
        before=tau*fraction
        time=tau*(fraction+1)
        E,n=exposure(before+u,arb(0),logs[:15])
        k1=before+q('.9999');k2=before+q('.0001')
        Z=(positive(k1)+positive(k2))/3
        def active(t):
            return arb(1) if t>0 else arb(0) if t<0 else arb(0).union(1)
        nf=(active(k1)+active(k2))/3
        xs=[v+arb(0,'.0015') for v in X];ys=[v+arb(0,'.003') for v in y]
        gs=couplings(xs,time,E,Z)
        xp,yp,fp=rates(gs,ys,n)
        gx=[p+nf*f for p,f in zip(xp,fp)]
        assert all(abs(dt*p)<q('.0015') for p in gx) and all(abs(dt*p)<q('.003') for p in yp)
        X=[v+dt*p for v,p in zip(X,gx)]
        y=[v+dt*p for v,p in zip(y,yp)]
    G=physical(X,tau,W,q('1/3'))
    ratio=rho
    # Common affine parts cancel exactly before interval arithmetic.
    D=X[0]-X[2]-2*((beta[0]-beta[2]+9*(eff[0]-eff[2]))*tau-9*(eff[0]-eff[2])*a+(hf[0]-hf[2])/3)/loop
    R=4*arb.pi()*(X[1]-ratio*X[0]-(1-ratio)*X[2]+18*(eff[1]-ratio*eff[0]-(1-ratio)*eff[2])*a/loop-2*(hf[1]-ratio*hf[0]-(1-ratio)*hf[2])/(3*loop))
    return R,D,u.exp(),4*arb.pi()*G[0]

faces=[('21.475','21.475','26.866','26.873'),('21.479','21.479','26.866','26.873'),
       ('21.475','21.479','26.866','26.866'),('21.475','21.479','26.873','26.873'),
       ('21.475','21.479','26.866','26.873')]
data=[run(*f) for f in faces]
assert data[0][0]<0 and data[1][0]>0 and data[2][1]>0 and data[3][1]<0
assert data[4][2]>q('456.75') and data[4][2]<q('459.79')
assert data[4][3]>q('38.3190') and data[4][3]<q('38.3265')
for row in data:print([(float(v.lower()),float(v.upper())) for v in row])
