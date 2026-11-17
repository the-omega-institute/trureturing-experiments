import sys
from fractions import Fraction as F
import sympy as sp

u,b,n,a = sp.symbols('u b n a')
THETA=sp.Rational(3,7); XI=sp.Rational(-35,29)
h=u*(u+1); ell=u*u+3*u-4
al=3*h-10; be=8*u*u+18*u-30
phi=al*a+be*n-23*h
delta=n*h+a*n*ell-4*a*h
A=(10-3*b)*u*u+(8-3*b)*u+10*b-32
B=(10-8*b)*u*u+(19-18*b)*u+30*b-36
C=(23*b-29)*h
Q=sp.Matrix([[1,-sp.Rational(1,2),-u-1],[-sp.Rational(1,2),-1,2*u+1],[-u-1,2*u+1,-u*u+u+1]])
T=sp.Matrix([[-al,(THETA*al-be-5*XI*ell)/2,(23+20*XI)*h/2],[(THETA*al-be-5*XI*ell)/2,THETA*be,-(23*THETA+5*XI)*h/2],[(23+20*XI)*h/2,-(23*THETA+5*XI)*h/2,0]])
vec=sp.Matrix([A,B,C])
rho=sp.expand((vec.T*Q.adjugate()*vec)[0]); sigma=sp.expand((vec.T*T.adjugate()*vec)[0])
tau=sp.expand((vec.T*((Q+T).adjugate()-Q.adjugate()-T.adjugate())*vec)[0])
R=sp.Poly(tau*tau-4*rho*sigma,u,b)
xa=sp.Matrix([-B*n-C,A*n,A])
qa=sp.Poly(sp.expand((xa.T*Q*xa)[0]),n); ta=sp.Poly(sp.expand((xa.T*T*xa)[0]),n)
q2,q1,q0=qa.all_coeffs(); t2,t1,t0=ta.all_coeffs()
nhat=sp.cancel(-sp.Rational(203,1)*(t2*q0-q2*t0)/A)
dhat=sp.cancel(-sp.Rational(203,1)*(q2*t1-t2*q1)/A)

def vals(ub, root):
    av=sp.N(A.subs({u:root,b:ub}),40); dv=sp.N(dhat.subs({u:root,b:ub}),40)
    nv=sp.N(nhat.subs({u:root,b:ub})/dhat.subs({u:root,b:ub}),40)
    aa=sp.N(-(B*nhat/dhat+C).subs({u:root,b:ub}),40)/av*av if False else sp.N((-(B*(nhat/dhat)-C*0 + C)/A).subs({u:root,b:ub}),40)
    # a=-(B*n+C)/A
    aa=sp.N((-(B*(nhat/dhat)+C)/A).subs({u:root,b:ub}),40)
    source=[nv+aa-3*root-1,nv+2*aa-4*root-3,1+2*root-aa,2+3*root-nv-aa]
    ph=sp.N(phi.subs({u:root,n:nv,a:aa}),40)
    drift=sp.N(120*root*(root+1)*XI/ph,40)
    g7=sp.N(210/nv+drift*((20*root**2+34*root-70)/(root*(root+1))-35/nv),40)
    mm=sp.N(210*(2+3*THETA)+30*XI*((4*(2*nv+3*aa)*(20*root**2+34*root-70)-105*aa*(root**2+3*root-4)-381*root*(root+1))/ph),40)
    mj=sp.N(drift+(4*nv-root+1)*g7,40)
    return nv,aa,source,mm,mj

def scan(bs):
    for bv in bs:
        poly=sp.Poly(R.as_expr().subs(b,sp.Rational(str(bv))),u)
        roots=sp.nroots(poly,maxsteps=200,n=35)
        out=[]
        for rr in roots:
            if abs(float(sp.im(rr)))<1e-18 and float(sp.re(rr))>0:
                ur=sp.re(rr)
                try:
                    nv,aa,ss,mm,mj=vals(sp.Rational(str(bv)),ur)
                    good=all(float(sp.re(x))>0 for x in ss)
                    if good:
                        out.append((float(ur),float(nv),float(aa),float(mm),float(mj)))
                except Exception:
                    pass
        print('b',bv,'worlds',len(out))
        for x in out: print(' ',tuple(round(y,10) for y in x))
        for name,idx in [('MM',3),('MJ',4)]:
            pairs=[]
            for i in range(len(out)):
                for j in range(i+1,len(out)):
                    if abs(out[i][idx]-out[j][idx])<1e-7: pairs.append((i+1,j+1))
            print(' ',name,'near_pairs',pairs)

if __name__=='__main__':
    scan([-100,-20,-5,-1,0,1,2,3,4,5,10,20,100])
