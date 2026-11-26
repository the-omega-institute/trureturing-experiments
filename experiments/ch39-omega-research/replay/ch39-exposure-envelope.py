"""Two-loop gauge enclosure for all threshold shapes with known maximum logarithmic widths."""
from pathlib import Path
from fractions import Fraction
from flint import arb
root=Path(__file__).parent
src=(root/'ch39-two-gauge-certificate.py').read_text().split('def run(al,')[0]
g={'__name__':'coefficients'}
exec(compile(src,'gauge_coefs','exec'),g)
q=g['q'];loop=g['loop'];beta=g['beta'];eff=g['eff'];hf=g['hf'];rho=g['rho']
pos=g['positive'];BD=g['BD'];BF=g['BF'];ref=g['gref'];gf=g['gf']
def mult(M,z):
    return [-2*sum(b*v for b,v in zip(row,z))/loop**2 for row in M]
def as_box(a,b):
    return q(a).union(q(b))
def positive_range(tau,rel,W,width,max_modes):
    if rel+width<0:return arb(0)
    # Every term vanishes at the left endpoint T-width and is convex.
    hi=(W*pos(1+rel/width)).upper()
    lo=pos(W+max_modes*rel).lower()
    return lo.union(hi)
def go(a0,a1,t0,t1,steps=2000):
    tau=as_box(t0,t1);a=as_box(a0,a1)
    W=9*(tau-a)
    assert W<52
    # The all-inertia width bound for W<=52; for smaller W verify tighter estimates.
    width=q('3.44')
    nsmax=38+270  # 2*135 Dirichlet first levels for R<1.62; W>54 at R=1.62.
    max_ferm=21
    X=g['G0'][:];y=g['y0'][:]
    dt=tau/steps
    MS=[None]*3;MF=[None]*3
    fixed_gs=as_box('.16','.58')
    rs=mult(BD,[fixed_gs-k for k in ref])
    rf=mult(BF,[fixed_gs-k for k in gf])
    for j in range(steps):
        rel=tau*as_box(Fraction(j-steps,steps),Fraction(j+1-steps,steps))
        time=tau*as_box(Fraction(j,steps),Fraction(j+1,steps))
        E=positive_range(tau,rel,W,width,nsmax)
        Z=positive_range(tau,rel,q('1/3'),q(1),q(max_ferm)/3)
        xs=[x+arb(0,'.002') for x in X];ys=[t+arb(0,'.006') for t in y]
        xx=[x+v*E+w*Z for x,v,w in zip(xs,rs,rf)]
        gc=g['couplings'](xx,time,E,Z)
        if not rel+width<0:
            assert all(fixed_gs.lower()<v and v<fixed_gs.upper() for v in gc)
            cs=mult(BD,[v-k for v,k in zip(gc,ref)])
            for i,v in enumerate(cs):MS[i]=v if MS[i] is None else MS[i].union(v)
        if not rel+1<0:
            cf=mult(BF,[v-k for v,k in zip(gc,gf)])
            for i,v in enumerate(cf):MF[i]=v if MF[i] is None else MF[i].union(v)
        xx=[x+(b or arb(0))*E+(c or arb(0))*Z for x,b,c in zip(xs,MS,MF)]
        gs=g['couplings'](xx,time,E,Z)
        xp,yp,_=g['rates'](gs,ys,0)
        assert all(abs(dt*v)<q('.002') for v in xp) and all(abs(dt*v)<q('.006') for v in yp)
        X=[x+dt*v for x,v in zip(X,xp)];y=[t+dt*v for t,v in zip(y,yp)]
    xr=[x+b*W+c/3 for x,b,c in zip(X,MS,MF)]
    G=g['physical'](xr,tau,W,q('1/3'))
    R=4*arb.pi()*(xr[1]-rho*xr[0]-(1-rho)*xr[2]+18*(eff[1]-rho*eff[0]-(1-rho)*eff[2])*a/loop-2*(hf[1]-rho*hf[0]-(1-rho)*hf[2])/(3*loop))
    D=xr[0]-xr[2]-2*((beta[0]-beta[2]+9*(eff[0]-eff[2]))*tau-9*(eff[0]-eff[2])*a+(hf[0]-hf[2])/3)/loop
    return R,D,4*arb.pi()*G[0],W

def certify(a,b,c,d):
    faces=((a,a,c,d),(b,b,c,d),(a,b,c,c),(a,b,d,d),(a,b,c,d))
    ans=[go(*face) for face in faces]
    for row in ans:print([(float(w.lower()),float(w.upper())) for w in row])
    assert ans[0][0]<0 and ans[1][0]>0 and ans[2][1]>0 and ans[3][1]<0
    return ans
if __name__=='__main__':
    import sys
    certify(*(sys.argv[1:] or ['21.43','21.52','26.837','26.899']))
