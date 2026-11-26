"""Clock phases at fixed r equal to the extreme arc ratio; certified continuum and finite-mesh estimates."""
from pathlib import Path
from fractions import Fraction
import runpy
from flint import arb, acb, acb_mat

root=Path(__file__).resolve().parent
flow=runpy.run_path(str(root/'ch39-spin-balanced-rge.py'))
q,D=flow['q'],flow['D']
data=flow['m']['data']; shell=data['shell']
r=data['theta'](804)/data['theta'](12)
pi=arb.pi()
cos_axis=(q(8741547)/12563572).sqrt()
sin_axis=(1-cos_axis**2).sqrt()
tF=data['m']['model']['geometry'](r)
d=data['d']
# 26 bounds every low endpoint Ritz energy (in fact <=16).
k=pi*(q(26)/d).sqrt()
eta=18*k*(k/2).tan()/26
# Removing the high sites above 959990 costs <= ||K_LH||^2/gap.
# Every row and column of this transport incidence block has at most 18 unit blocks.
high_error=q(18)**2/(q(959990)-26)
def low_matrix(p,second=None,axis=None):
    sine,cosine=(sin_axis,cos_axis) if axis is None else axis.sin_cos()
    s,c=p.sin_cos()
    U=acb_mat([[acb(c,s),0],[0,acb(c,-s)]])
    v,w=(-p if second is None else second).sin_cos()
    V=acb_mat([[acb(w,cosine*v),sine*v],[-sine*v,acb(w,-cosine*v)]])
    H=acb_mat(18,18)
    for i in range(9):
        internal=sum(i in (a,b) for a,b,_,_ in shell['edges'] if b<9)
        for h in range(2):H[2*i+h,2*i+h]=internal
    for i,j,ds,z in shell['edges']:
        if j>=9:continue
        A=U if ds==724 else V if ds==804 else acb_mat([[1,0],[0,1]])
        for a in range(2):
            for b in range(2):
                val=-z*A[a,b]
                H[2*i+a,2*j+b]=val
                H[2*j+b,2*i+a]=val.conjugate()
    return H

def inertia(H,x):
    K=acb_mat(H)
    for i in range(18):K[i,i]-=x
    ids=list(range(18))
    n=0
    while ids:
        j=max(ids,key=lambda j:float(abs(K[j,j].real).mid()))
        ids.remove(j)
        pivot=K[j,j].real
        assert pivot<0 or pivot>0
        n+=int(pivot<0)
        for h,i in enumerate(ids):
            t=K[i,j]/pivot
            for k in ids[h:]:
                v=K[i,k]-t*K[j,k]
                K[i,k]=v.real if i==k else v
                if i!=k:K[k,i]=v.conjugate()
    return n

def spectrum(p,width,accurate=False,second=None,axis=None):
    H=low_matrix(p,second,axis)
    # Build an ordered enclosure recursively, sharing inertia counts between adjacent eigenvalues.
    ans=[None]*18
    def isolate(a,b,na,nb,depth):
        if na==nb:return
        if depth==0:
            for j in range(na,nb):ans[j]=q(a).union(q(b))
            return
        midpoint=(a+b)/2
        count=inertia(H,q(midpoint))
        isolate(a,midpoint,na,count,depth-1)
        isolate(midpoint,b,count,nb,depth-1)
    a,b=Fraction(1,4),Fraction(65,4)
    assert inertia(H,q(a))==0 and inertia(H,q(b))==18
    isolate(a,b,0,18,15 if accurate else 7)
    # The maximal number of clocked low-low edges incident at a low site is four.
    assert max(sum(i in (a,b) and b<9 and ds in (724,804) for a,b,ds,_ in shell['edges']) for i in range(9))==4
    ys=[]
    for interval in ans:
        lo=(interval.lower()-4*width-high_error)/(1+eta)
        lo=lo/(1+data['beta']*(lo+10))
        hi=interval.upper()+4*width
        assert lo>0 and hi<26
        ys.append(lo.lower().union(hi.upper()))
    return ys

def run(tau,xs,steps=1024,scale=r,return_target=False):
    g=flow['g']; m=flow['m']
    tF=data['m']['model']['geometry'](scale)
    logs=[scale.log()-x.log()/2 for x in xs]
    w=sum(logs)
    X=g['G0'][:];Y=g['y0'][:]
    padX=q('.001')*1024/steps;padY=q('.006')*1024/steps
    dt=tau/steps
    for j in range(steps):
        rel=tau*q(Fraction(j-steps,steps)).union(q(Fraction(j+1-steps,steps)))
        t=tau*q(Fraction(j,steps)).union(q(Fraction(j+1,steps)))
        zs=[rel+z for z in logs]
        el=sum(g['positive'](z) for z in zs);nl=sum(m['active'](z) for z in zs)
        ef=sum(g['positive'](rel+z) for z in tF)/3; nf=sum(m['active'](rel+z) for z in tF)/3
        tube=[v+arb(0,padX) for v in X];ys=[v+arb(0,padY) for v in Y]
        G=[v-2*(b*t+e*el+f*ef)/D for v,b,e,f in zip(tube,flow['beta'],flow['eff'],g['hf'])]
        assert all(v>1/(4*pi) for v in G) and all(y.exp()<4*pi for y in ys)
        gs=[1/v for v in G];shift=[v-z for v,z in zip(gs,g['gref'])]
        xp,yp,fp=g['rates'](gs,ys,q(0))
        xp=[v-2*nl*m['dot'](row,shift)/D**2+nf*f for v,row,f in zip(xp,flow['B'],fp)]
        assert all(abs(dt*v)<padX for v in xp) and all(abs(dt*v)<padY for v in yp)
        X=[v+dt*f for v,f in zip(X,xp)];Y=[v+dt*f for v,f in zip(Y,yp)]
    target=[D*x/2-f/3 for x,f in zip(X,g['hf'])]
    inv=4*pi*(X[0]-2*(flow['beta'][0]*tau+flow['eff'][0]*w+g['hf'][0]/3)/D)
    time=m['dot'](flow['pt'],target)
    return (time if return_target else time-tau),m['dot'](flow['pw'],target)-w,w,inv

if __name__=='__main__':
    t0,t1=q('26.94'),q('27.06')
    T=t0.union(t1)
    a,b=3*pi/8,pi
    for phase,direction in ((a,1),(b,-1)):
        xs=spectrum(phase,q(0),True)
        pred,wres,w,_=run(T,xs,4096)
        print('endpoint log exposure',w.lower().str(12),w.upper().str(12),
              'signed exposure residual', (direction*wres).lower().str(12),flush=True)
        assert direction*wres>0
    pieces=16
    margins=[]
    for j in range(pieces):
        width=(b-a)/(2*pieces)
        xs=spectrum(a+(2*j+1)*width,width)
        left=run(t0,xs)[0];right=run(t1,xs)[0]
        margins.append(min(left.lower(),(-right).lower()))
    print('opposite time face margin',min(margins).str(12))
    assert min(margins)>0
