"""Two-loop gauge and one-loop light Yukawa crossing for spin-balanced isospectral towers."""
from pathlib import Path
from fractions import Fraction
import runpy
from flint import arb

root = Path(__file__).resolve().parent
m = runpy.run_path(str(root/'ch39-right-scalar-rge.py'))
g = m['g']
q,D = g['q'],g['loop']
L,R = g['BD'],m['BR']
def dirac(i,j):
    return sum(T[i]*(4*C[j]+q('20/3')*(0,2,3)[i]*int(i==j))
               for T,C in zip(g['indices'],g['casimir']))
B = [[L[i][j]+3*R[i][j]+dirac(i,j) for j in range(3)] for i in range(3)]
beta = g['beta']
d = [2+2*t for t in g['d']]
eff = [di+m['dot'](row,g['gref'])/D for di,row in zip(d,B)]
pt,pw = m['project'](eff,beta),m['project'](beta,eff)
low = m['data']['xs']
shape = [x.log()/2 for x in low]
def projected(tau,w,steps=4096):
    logr = (w+sum(shape))/18
    r = logr.exp()
    assert r**2<m['data']['edge19'] and all(x<r**2 for x in low)
    # Work with the total W as a coordinate to avoid counting low-spectrum uncertainties twice.
    exposures = [logr-a for a in shape]
    assert all(z<tau for z in exposures)
    tF = m['data']['m']['model']['geometry'](r)
    X = g['G0'][:]
    Y = g['y0'][:]
    dt = tau/steps
    for j in range(steps):
        rel = tau*q(Fraction(j-steps,steps)).union(q(Fraction(j+1-steps,steps)))
        t = tau*q(Fraction(j,steps)).union(q(Fraction(j+1,steps)))
        active = [rel+li for li in exposures]
        el = sum(g['positive'](z) for z in active)
        nl = sum(m['active'](z) for z in active)
        ef = sum(g['positive'](rel+v) for v in tF)/3
        nf = sum(m['active'](rel+v) for v in tF)/3
        tube = [v+arb(0,'.001') for v in X]
        ys = [v+arb(0,'.006') for v in Y]
        G = [v-2*(b*t+e*el+f*ef)/D for v,b,e,f in zip(tube,beta,eff,g['hf'])]
        assert all(v>1/(4*arb.pi()) for v in G) and all(y.exp()<4*arb.pi() for y in ys)
        gs = [1/v for v in G]
        xp,yp,fp = g['rates'](gs,ys,q(0))
        shift = [v-r0 for v,r0 in zip(gs,g['gref'])]
        xp = [v-2*nl*m['dot'](row,shift)/D**2+nf*f for v,row,f in zip(xp,B,fp)]
        assert all(abs(dt*v)<q('.001') for v in xp) and all(abs(dt*v)<q('.006') for v in yp)
        X = [v+dt*k for v,k in zip(X,xp)]
        Y = [v+dt*k for v,k in zip(Y,yp)]
    target = [D*x/2-f/3 for x,f in zip(X,g['hf'])]
    alpha = 4*arb.pi()*(X[0]-2*(beta[0]*tau+eff[0]*w+g['hf'][0]/3)/D)
    return m['dot'](pt,target),m['dot'](pw,target),r,alpha
def run(tau,w,steps=4096):
    tstar,wstar,r,alpha = projected(tau,w,steps)
    return tstar-tau,wstar-w,r,alpha
if __name__=='__main__':
    import sys
    t0,t1,w0,w1=map(q,sys.argv[1:] or ('26.99','27.013','24.925','25.07'))
    T,W=t0.union(t1),w0.union(w1)
    faces=[run(t0,W),run(t1,W),run(T,w0),run(T,w1),run(T,W)]
    for i,row in enumerate(faces):
        print('face',i,[(v.lower().str(10),v.upper().str(10)) for v in row])
    assert faces[0][0]>0 and faces[1][0]<0 and faces[2][1]>0 and faces[3][1]<0
