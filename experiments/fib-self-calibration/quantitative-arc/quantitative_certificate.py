from fractions import Fraction as Q
from math import comb
class I:
    __slots__=('lo','hi')
    def __init__(self,a=0,b=None):
        if isinstance(a,I): self.lo,self.hi=a.lo,a.hi; return
        self.lo=Q(a); self.hi=self.lo if b is None else Q(b); assert self.lo<=self.hi
    def __add__(self,y): y=I(y); return I(self.lo+y.lo,self.hi+y.hi)
    __radd__=__add__
    def __neg__(self): return I(-self.hi,-self.lo)
    def __sub__(self,y): return self+-I(y)
    def __rsub__(self,y): return I(y)+-self
    def __mul__(self,y):
        y=I(y); z=[self.lo*y.lo,self.lo*y.hi,self.hi*y.lo,self.hi*y.hi]; return I(min(z),max(z))
    __rmul__=__mul__
    def __truediv__(self,y):
        y=I(y); assert y.hi<0 or y.lo>0, (y.lo,y.hi); return self*I(1/y.hi,1/y.lo)
    def __rtruediv__(self,y): return I(y)/self
    def __pow__(self,n):
        z=I(1)
        for _ in range(n): z=z*self
        return z
    def intersect(self,y):
        y=I(y); a=max(self.lo,y.lo); b=min(self.hi,y.hi); return None if a>b else I(a,b)
    def nz(self): return self.hi<0 or self.lo>0
    def width(self): return self.hi-self.lo

def Hn(c,x):
    z=I()
    for a in c: z=z*x+a
    return z

def polyval(a,t,b): return Hn([Hn(row[::-1],b) for row in a[::-1]],t)

def j_floor():
    b=I(0,Q(1,1024)); z=I(0,Q(1,1024)); T=1+b; U=1+2*b; q=-6*(8*b-1)*(5*b+2)/(5*T*U); out={}
    for mode in ('low','high'):
        if mode=='low':
            u=I(0,Q(1,2048)); n=z*(1+u); den=z*(1+6/T+1/U)-4-4/(1+u); w=q/den; e=z*w
            vals=[-24-e,-w,(Q(27,5)-e)/T,-w/(1+u)]
        else:
            a=I(0,Q(1,2*1024**2)); n=z+1-a; den=z*(1+6/T+1/U)-4-4*z/n; w=q/den; e=z*w
            vals=[-24-e,-w,(Q(27,5)-e)/T,-z*w/n]
        g5=sum((a*c for a,c in zip([-4,15,-20,10],vals)),I()); g6=sum((a*c for a,c in zip([-10,36,-45,20],vals)),I())
        p=24+e+(U+n)*g5; v=24+e+(2*U+n)*g6; gap=v-Q(37,10)*p-Q(213,5)
        assert den.nz() and gap.hi < -20; out[mode]=(den,gap)
    return out

def source_cover():
    delta=Q(1,10**6); stack=[(I(0,1),I(0,30),I(Q(1,1024),40),0)]; visits=discarded=retained=splits=depthmax=0; reasons={}
    def toss(k): reasons[k]=reasons.get(k,0)+1
    while stack:
        rr,bb,jj,depth=stack.pop(); visits+=1; depthmax=max(depthmax,depth)
        if rr.lo>=delta and rr.hi<=1-delta and bb.lo>=delta and jj.lo-bb.hi>=delta: retained+=1; continue
        J=jj.intersect(bb+I(0,10))
        if J is None: discarded+=1; toss('c_bound'); continue
        qlo=min(rr.lo*(1-rr.lo),rr.hi*(1-rr.hi)); qhi=Q(1,4) if rr.lo<=Q(1,2)<=rr.hi else max(rr.lo*(1-rr.lo),rr.hi*(1-rr.hi)); qr=I(qlo,qhi)
        if bb.lo>0:
            J=J.intersect(bb+qr/bb)
            if J is None: discarded+=1; toss('rank_contract'); continue
        if (bb*(J-bb)-qr).nz(): discarded+=1; toss('rank_residual'); continue
        T=1+bb; U=1+2*bb; n=J+rr; la=1+6/T+1/U-4/J-4/n; db=6*(8*bb-1)*(5*bb+2)/(5*T*U)
        p0=24+(U+n)*(96-108/T); ap=1+(U+n)*(4-15/J+20/T-10/n); v0=24+(2*U+n)*(240-243/T); bv=1+(2*U+n)*(10-36/J+45/T-20/n)
        E=I(-655000,655000); failed=False
        for aa,hh in ((la,-db),(ap,I(10,34)-p0),(bv,I(0,200)-v0)):
            if aa.nz():
                E=E.intersect(hh/aa)
                if E is None: failed=True; break
        if failed: discarded+=1; toss('e_contract'); continue
        if any((aa*E-hh).nz() for aa,hh in ((la,-db),(ap,I(10,34)-p0),(bv,I(0,200)-v0))): discarded+=1; toss('response_residual'); continue
        gap=v0-Q(37,10)*p0-Q(213,5)+(bv-Q(37,10)*ap)*E
        if gap.hi<0: discarded+=1; toss('barrier'); continue
        widths=[rr.width(),bb.width()/3,jj.width()/4]; ix=widths.index(max(widths)); boxes=[rr,bb,jj]; ii=boxes[ix]; mid=(ii.lo+ii.hi)/2; assert ii.lo<mid<ii.hi
        lo=boxes.copy(); hi=boxes.copy(); lo[ix]=I(ii.lo,mid); hi[ix]=I(mid,ii.hi); stack.append((*hi,depth+1)); stack.append((*lo,depth+1)); splits+=1
    assert visits==2*(discarded+retained)-1 and splits==discarded+retained-1
    return {'visits':visits,'discarded':discarded,'retained':retained,'pending':len(stack),'max_depth':depthmax,'splits':splits,'reasons':reasons}

blocks=[
(10,20,[-286517520774,-1458798946604,-4484619212987,-39378923335258,-95032360161001,-610511414047453,-3787735814587272,5327878585599447,-134291913516490829,725907457857128599,3554914500654767298,210788533712955551542],[-21554952211,-167568989307,-719043316224,-5332188024970,-17244126711789,-76110792776341,-367687512419822,1518778661443952,1190131346284344,224217790257326402,4016662099446683744,4532053336144848092]),
(20,30,[-64760914995336,-141100470063114,-72520370219354,-145269729373691,-728244036456201,-2245702071706101,-7245822024720854,-23692327394115831,-73163432013250704,-197029754390717363,-528694944588356620,-584207156211347097,4492542348396760978,219625886685223990790],[-1723603413335,-1880359411512,11905805971444,49433821136737,124147171707093,215207360044166,34937022519674,-1562471228161450,-7777580962203138,-22939857428839428,-39711075675135058,193301006165304340,4916488201945917316,13475282226511291486]),
(30,34,[-55541680313275,-97404534628203,106298090982829,184220076471718,-270290783930968,-478620521827863,-157686308174145,-308182157741107,-1102948577953311,-2087793427689654,-3803786492746523,-7536086716159247,-15362660199585885,-32003595231890594,-70148176855102684,-175690543223014516,-538164996818189229,-1728970444713086177,-2173832085911805516,221845194568232864125],[450135834043,791432409571,-853929166848,-1481671306073,2207327774598,3924165644115,1387884608488,2738482983299,9556425419917,18621330504658,34279898879426,63331218908555,127163137096159,401052794968752,1494889997020038,2689161815711592,-12134528603054534,-66100473199783546,1928352859377411297,20494322221751046893])]
Darr=[[0,-20,-40,0,5],[-20,-30,30,30,5],[20,80,65,-35,-10],[0,-40,-100,-50,-5],[0,0,20,40,5],[0,0,0,0,0],[0,0,0,0,0]]
Varr=[[0,-1032,-7458,-13293,-1125],[-1032,-7320,-11472,9750,10650],[750,6891,26640,25266,-9750],[270,1695,-5574,-33168,-17880],[0,-810,-8295,-5595,13050],[0,0,810,8505,5310],[0,0,0,-270,-2655]]
rho=Q(1,100000); arc={}
for lo,hi,tc,bc in blocks:
    vlo=Q(10**9); vhi=-vlo; dmin=Q(10**9); dmax=-vlo
    for k in range(1024):
        u0=Q(-1)+Q(2*k,1024); u1=Q(-1)+Q(2*(k+1),1024)
        def centre(cs,x): return Hn([Q(a,10**20) for a in cs],x)
        tb=centre(tc,I(u0,u1))+I(-rho,rho); bb=centre(bc,I(u0,u1))+I(-rho,rho); di=polyval(Darr,tb,bb); vi=polyval(Varr,tb,bb)/di
        assert di.lo>0; dmin=min(dmin,di.lo); dmax=max(dmax,di.hi); vlo=min(vlo,vi.lo); vhi=max(vhi,vi.hi)
    arc[(lo,hi)]=(dmin,dmax,vlo,vhi); simple={(10,20):(Q(81),Q(118)),(20,30):(Q(116),Q(156)),(30,34):(Q(153),Q(171))}[(lo,hi)]; assert simple[0]<vlo and vhi<simple[1]

import sympy as s
t,b,r,ss,b1,b2=s.symbols('t b r ss b1 b2'); j=t+b*(1-t*t); n=j+t*b; T=1+b; U=1+2*b; D=5*((b*b+8*b+4)*j*n-2*T*U*(j+n)); A=(8*b-1)*(2+5*b)
C0=-24*D+3*A*j*n; C1=3*A*n; C2=3*(49*b+34)*j*n-54*U*(j+n); C3=3*A*j; B5=-4*C0+15*C1-20*C2+10*C3; B6=-10*C0+36*C1-45*C2+20*C3; B7=-20*C0+70*C1-84*C2+35*C3
Np=s.expand(-C0+(U+n)*B5); Nv=s.expand(-C0+(2*U+n)*B6); Nm=s.expand(-C0+(3*U+2*n)*B7); Snum=s.Poly(-(Np-24*D),t,b); Cnum=s.Poly(-(5*Nm-32*Nv+1128*D),t,b); den=s.Poly(3*(8*b-1),t,b); Sq,Sr=s.div(Snum,den); Cq,Cr=s.div(Cnum,den); assert Sr.is_zero and Cr.is_zero
H=14*b**6+(95-73*r)*b**5+(47*r**2-39*r+30)*b**4+r*(147*r**2-232*r+105)*b**3+2*r*(1-r)*(44*r**2-55*r+15)*b**2-64*r**3*(1-r)**2*b-32*r**3*(1-r)**3
Cstar=-1055*b**6+(5815*r-7529)*b**5+(-5140*r**2+5260*r-3102)*b**4+(-8795*r**3+13798*r**2-6351*r)*b**3+(7015*r**4-16924*r**3+13011*r**2-3102*r)*b**2+(1920*r**5-1610*r**4-2540*r**3+2230*r**2)*b-960*r**6+2880*r**5-2880*r**4+960*r**3
assert s.simplify(b**3*Sq.as_expr().subs({t:r/b})-H)==0 and s.simplify(b**3*Cq.as_expr().subs({t:r/b})-Cstar)==0
F=s.expand(Cstar.subs({r:ss,b:b2})*H.subs({r:r,b:b1})-Cstar.subs({r:r,b:b1})*H.subs({r:ss,b:b2})); F0=21*(16*r-9)*(64*r*r-72*r-1)*(16777216*r**6-56623104*r**5+93622272*r**4-91210752*r**3+37316032*r**2+622224*r-230227)/s.Integer(68719476736); assert s.simplify(F.subs({ss:s.Rational(9,8)-r,b1:s.Rational(1,8),b2:s.Rational(1,8)})-F0)==0

def bern(poly,var,a,bb):
    x=s.symbols('x'); q=s.Poly(s.expand(poly.subs(var,s.Rational(a.numerator,a.denominator)+(s.Rational(bb.numerator,bb.denominator)-s.Rational(a.numerator,a.denominator))*x)),x); n=q.degree(); c=[Q(int(q.nth(i).p),int(q.nth(i).q)) for i in range(n+1)]; return [sum(c[i]*Q(comb(k,i),comb(n,i)) for i in range(k+1)) for k in range(n+1)]
fb=bern(F0,r,Q(267,1000),Q(277,1000)); fmin=min(fb); assert fmin==Q(14656469065762145135535204561,512000000000000000000000000000) and fmin>Q(1,35)
def qsym(x):
    if isinstance(x,Q): return x
    if hasattr(x,'p') and hasattr(x,'q'): return Q(int(x.p),int(x.q))
    return Q(x)
def ieval(poly,vars,box):
    if not vars: return I(qsym(poly))
    P=s.Poly(s.expand(poly),vars[0]);
    if P.is_zero: return I()
    v=vars[0]
    cs=[ieval(P.nth(i),vars[1:],box) for i in range(P.degree(v),-1,-1)]
    return Hn(cs,box[v])

box={r:I(Q(267,1000),Q(277,1000)),ss:I(Q(9,8)-Q(277,1000)-Q(1,100000),Q(9,8)-Q(267,1000)+Q(1,100000)),b1:I(Q(1,8)-Q(1,100000),Q(1,8)+Q(1,100000)),b2:I(Q(1,8)-Q(1,100000),Q(1,8)+Q(1,100000))}
def two(poly,rr,bb):
    return ieval(poly,(r,b),{r:rr,b:bb})
H1=two(H,box[r],box[b1]); H2=two(H,box[ss],box[b2]); C1i=two(Cstar,box[r],box[b1]); C2i=two(Cstar,box[ss],box[b2])
mirror_derivatives={
 'b1':C2i*two(s.diff(H,b),box[r],box[b1])-two(s.diff(Cstar,b),box[r],box[b1])*H2,
 'b2':two(s.diff(Cstar,b),box[ss],box[b2])*H1-C1i*two(s.diff(H,b),box[ss],box[b2]),
 's':two(s.diff(Cstar,r),box[ss],box[b2])*H1-C1i*two(s.diff(H,r),box[ss],box[b2])}
for var,cap in (('b1',20),('b2',25),('s',65)):
    q=mirror_derivatives[var]; assert q.lo>-cap and q.hi<cap, (var,str(q.lo),str(q.hi))
assert fmin-Q(110,100000)>Q(1,40)

def bern1(cs,a,bb):
    x=s.symbols('x'); n=len(cs)-1; poly=sum(s.Rational(c.numerator,c.denominator)*(s.Rational(a.numerator,a.denominator)+(s.Rational(bb.numerator,bb.denominator)-s.Rational(a.numerator,a.denominator))*x)**i for i,c in enumerate(cs)); q=s.Poly(s.expand(poly),x); co=[Q(int(q.nth(i).p),int(q.nth(i).q)) for i in range(n+1)]; return [sum(co[i]*Q(comb(k,i),comb(n,i)) for i in range(k+1)) for k in range(n+1)]
P=[Q(-225625),Q(-21500400),Q(41464676),Q(34142752),Q(-53190656),Q(-13631488),Q(16777216)]; BB=[Q(-900),Q(31),Q(716)]; R=[Q(-57875),Q(-12320),Q(119952),Q(-8192)]; lo,hi=Q(47,24),Q(2); PB,BBB,RB=bern1(P,lo,hi),bern1(BB,lo,hi),bern1(R,lo,hi)
assert min(PB)==Q(842710954799,5832) and max(PB)==Q(182258007) and min(BBB)==Q(274553,144) and max(BBB)==Q(2026) and min(RB)==Q(34181767,108) and max(RB)==Q(331757)
QL=[48*R[i]-61*125*(BB[i] if i<len(BB) else 0) for i in range(4)]; QH=[65*125*(BB[i] if i<len(BB) else 0)-48*R[i] for i in range(4)]; assert min(bern1(QL,lo,hi))==Q(476086) and min(bern1(QH,lo,hi))==Q(43110037,144)
assert Q(57,100)**2<Q(21,64) and Q(67,192)<Q(3,5)**2
Plo,Phi=min(PB),max(PB); Blo,Bhi=min(BBB),max(BBB); Rlo,Rhi=min(RB),max(RB); beta_abs=Q(175)*I(lo,hi)*I(Blo,Bhi)*I(Q(57,100),Q(3,5))/I(Plo,Phi); mu=Q(336,5)*I(lo,hi)*I(Rlo,Rhi)*I(Q(57,100),Q(3,5))/I(Plo,Phi)
assert Q(1,600)<beta_abs.lo and beta_abs.hi<Q(1,250) and Q(1,10)<mu.lo and mu.hi<Q(1,4)

# Full-parent analytic derivatives, reduced in the quadratic mirror algebra.
N,zeta=s.symbols('N zeta')
rr=s.Rational(9,16)-zeta/2; sr=s.Rational(9,16)+zeta/2
mirror_relation=zeta**2-(s.Rational(85,32)-N)/2
Pn=sum(s.Rational(c.numerator,c.denominator)*N**i for i,c in enumerate(P))
Bn=sum(s.Rational(c.numerator,c.denominator)*N**i for i,c in enumerate(BB))
Rn=sum(s.Rational(c.numerator,c.denominator)*N**i for i,c in enumerate(R))
def reduction_zero(expr):
    num=s.cancel(expr).as_numer_denom()[0]
    assert s.Poly(num,zeta,domain=s.QQ.frac_field(N)).rem(s.Poly(mirror_relation,zeta,domain=s.QQ.frac_field(N))).is_zero

def jet(x):
    jx=N-x
    lam=s.Rational(107,15)-4/jx-4/N
    be=-s.Rational(25,448)*lam
    gamma5=(-256+85/jx+190/N)/21
    gamma6=(-185+48/jx+160/N)/7
    gamma7=-48+10/jx+45/N
    alpha=1+(N+s.Rational(5,4))*gamma5
    eta=-s.Rational(27,5)+(s.Rational(15,4)+2*N)*gamma7-s.Rational(32,5)*(s.Rational(5,2)+N)*gamma6-s.Rational(768,5)*be
    return be,alpha,eta
ber,ar,etar=jet(rr); bes,ass,etas=jet(sr)
reduction_zero(ar*ass-Pn/(1764*N**2*(64*N**2-64*N-1)))
reduction_zero(bes/ass-ber/ar+175*N*Bn*zeta/Pn)
reduction_zero(etas/ass-etar/ar-s.Rational(336,5)*N*Rn*zeta/Pn)
assert Q(21,64)>Q(571,1000)**2 and Q(67,192)<Q(591,1000)**2

# Zero-safe parent/comparator identities with independent j,n,e.
jj,nn,ee=s.symbols('jj nn ee')
TT=1+b; UU=1+2*b
lam=1+6/TT+1/UU-4/jj-4/nn
db=6*(8*b-1)*(5*b+2)/(5*TT*UU)
gs=[-24-ee,-ee/jj,(s.Rational(27,5)-ee)/TT,-ee/nn]
def row(cs): return sum(x*y for x,y in zip(cs,gs))
g5=row([-4,15,-20,10]); g6=row([-10,36,-45,20]); g7=row([-20,70,-84,35])
pp=24+ee+(UU+nn)*g5; vv=24+ee+(2*UU+nn)*g6; mm=24+ee+(3*UU+2*nn)*g7
parent=ee*lam+db
q0=s.factor((5*mm-32*vv+1128).subs(ee,0))
KC=s.cancel((5*b+2)*jj*nn*q0/db)
KS=10*UU*jj*nn*(UU+nn)
SG=s.cancel((5*b+2)*jj*nn*s.diff(pp,ee)-KS*lam)
CG=s.cancel((5*b+2)*jj*nn*s.diff(5*mm-32*vv,ee)-KC*lam)
assert s.cancel((5*b+2)*jj*nn*(pp-24)-ee*SG-KS*parent)==0
assert s.cancel((5*b+2)*jj*nn*(5*mm-32*vv+1128)-ee*CG-KC*parent)==0
assert s.cancel(SG.subs({jj:j,nn:n})-Sq.as_expr())==0
assert s.cancel(CG.subs({jj:j,nn:n})-Cq.as_expr())==0
assert s.cancel(KC--10*UU*jj*nn*(170*b+30*nn+89))==0

# Every coefficient used for the actual-arc enclosure is reconstructed.
for expr,arr in ((D,Darr),(Nv,Varr)):
    poly=s.Poly(s.expand(expr),t,b)
    assert [[int(poly.coeff_monomial(t**i*b**k)) for k in range(5)] for i in range(7)]==arr

# Exact characteristic-zero h divisibility; no sampled evaluation.
h,sigma=s.symbols('h sigma')
def shifted(poly): return poly.subs({t:t+h*sigma,b:b+h},simultaneous=True)
quotients={}
for name,f,g in [('p',Np,D),('v',Nv,D),('comparator',Cq.as_expr(),Sq.as_expr())]:
    poly=s.Poly(s.expand(shifted(f)*g-f*shifted(g)),h,t,b,sigma)
    assert poly.eval(h,0).is_zero
    quotient,remainder=s.div(poly,s.Poly(h,h,t,b,sigma))
    assert remainder.is_zero
    quotients[name]={'monomials':len(quotient.terms()),'total_degree':quotient.total_degree(),'remainder':0}
assert quotients=={'p':{'monomials':1518,'total_degree':23,'remainder':0},'v':{'monomials':1518,'total_degree':23,'remainder':0},'comparator':{'monomials':853,'total_degree':21,'remainder':0}}

# Reuse of the published scalar barrier, with exact rational MVT arithmetic.
def derivative_array(a,axis):
    out=[[0]*5 for _ in range(7)]
    for i in range(7):
        for k in range(5):
            if axis==0 and i: out[i-1][k]=i*a[i][k]
            if axis==1 and k: out[i][k-1]=k*a[i][k]
    return out
A0=[[10*Varr[i][k]-426*Darr[i][k] for k in range(5)] for i in range(7)]
A1=[[-37*x for x in row] for row in Darr]
def affine_eval(a0,a1,tb,bb,pb):
    return polyval([[I(x)+pb*y for x,y in zip(row0,row1)] for row0,row1 in zip(a0,a1)],tb,bb)
Dt0,Dt1=derivative_array(A0,0),derivative_array(A1,0)
Db0,Db1=derivative_array(A0,1),derivative_array(A1,1)
barriers={}
for low,high,tc,bc in blocks:
    minimum=None; min_index=None
    for k in range(1024):
        uk=Q(-1)+Q(2*k+1,1024); ui=I(uk-Q(1,1024),uk+Q(1,1024))
        ct=Hn([Q(x,10**20) for x in tc],I(uk)); cb=Hn([Q(x,10**20) for x in bc],I(uk))
        xt=Hn([Q(x,10**20) for x in tc],ui)+I(-rho,rho); xb=Hn([Q(x,10**20) for x in bc],ui)+I(-rho,rho)
        pc=Q(low+high,2)+Q(high-low,2)*uk; pi=I(Q(low+high,2))+Q(high-low,2)*ui
        rc=affine_eval(A0,A1,ct,cb,I(pc))
        ht=affine_eval(Dt0,Dt1,xt,xb,pi); hb=affine_eval(Db0,Db1,xt,xb,pi); hp=polyval(A1,xt,xb)
        displacements=[max(ct.hi-xt.lo,xt.hi-ct.lo),max(cb.hi-xb.lo,xb.hi-cb.lo),max(pc-pi.lo,pi.hi-pc)]
        bound=rc.lo-sum(max(abs(q.lo),abs(q.hi))*delta for q,delta in zip((ht,hb,hp),displacements))
        assert bound>0
        if minimum is None or bound<minimum: minimum,min_index=bound,k
    barriers[str((low,high))]={'cells':1024,'min_index':min_index,'lower':str(Q((minimum*10**12).__floor__(),10**12))}

# The published p=24 endpoint enclosure and actual mirror correspondence.
_,_,tc,bc=blocks[1]
pt=Hn([Q(x,10**20) for x in tc],I(Q(-1,5)))+I(-rho,rho)
pb=Hn([Q(x,10**20) for x in bc],I(Q(-1,5)))+I(-rho,rho)
p24=polyval(Varr,pt,pb)/polyval(Darr,pt,pb)
assert p24.lo>131 and p24.hi<132

import json
from pathlib import Path

def outward(x,upper=False,scale=10**6):
    y=x*scale
    return str(Q(-((-y.numerator)//y.denominator) if upper else y.numerator//y.denominator,scale))
def interval_json(x): return [str(x.lo),str(x.hi)]
cover_result=source_cover()
assert cover_result=={'visits':169103,'discarded':52765,'retained':31787,'pending':0,'max_depth':49,'splits':84551,'reasons':{'barrier':1224,'rank_residual':14136,'rank_contract':4901,'c_bound':4724,'e_contract':27652,'response_residual':128}}

# Optional exact comparison with the mathematical polynomial expressions.
import argparse,re,hashlib
from sympy.parsing.sympy_parser import parse_expr,standard_transformations,implicit_multiplication_application,convert_xor
parser=argparse.ArgumentParser()
parser.add_argument('--theory',type=Path)
args=parser.parse_args()
theory_correspondence=None
if args.theory:
    raw=args.theory.read_bytes()
    assert hashlib.sha256(raw[:48060]).hexdigest()=='fca9eccbaa14fa415fbfa57ef718f6f0d5b52d793d287f1251fd9acee9d2a93c'
    text=raw[48060:].decode()
    def mathematical_poly(fragment):
        fragment=re.sub(r'\\frac\{([^{}]*)\}\{([^{}]*)\}',r'((\1)/(\2))',fragment)
        for macro in ('\\bigl','\\bigr','\\quad','\\qquad'):
            fragment=fragment.replace(macro,'')
        fragment=fragment.replace('\\\\','').replace('&','').replace('{','').replace('}','').replace('\n','').strip().rstrip(',.;')
        return parse_expr(fragment,local_dict={'r':r,'b':b,'n':N},transformations=standard_transformations+(implicit_multiplication_application,convert_xor))
    ht=text.split('\\mathscr H(r,b)={}&',1)[1].split('C_*(r,b)={}&',1)[0]
    ct=text.split('C_*(r,b)={}&',1)[1].split('\\end{aligned}',1)[0]
    ft=text.split('F_0(r)=',1)[1].split('\\end{aligned}',1)[0]
    pt=text.split('P(n)={}&',1)[1].split('B(n)={}&',1)[0]
    bt=text.split('B(n)={}&',1)[1].split('R(n)={}&',1)[0]
    rt=text.split('R(n)={}&',1)[1].split('\\end{aligned}',1)[0]
    for fragment,poly in ((ht,H),(ct,Cstar),(ft,F0),(pt,Pn),(bt,Bn),(rt,Rn)):
        assert s.expand(mathematical_poly(fragment)-poly)==0
    for value in ('169103','84551','31787','4724','4901','14136','27652','128','1224','1518','853','1923','70000'):
        assert value in text
    for var,low,high in [('b1',-11,19),('b2',-24,19),('s',-62,58)]:
        assert low<mirror_derivatives[var].lo and mirror_derivatives[var].hi<high
    assert text.endswith('## 追加锚（本行以下为增补区）\n')
    assert not any(x in text for x in ('~~~','```','\\(','\\[','sshx','attempt-','/Users/'))
    assert len(raw.splitlines())<=5000
    theory_correspondence={'polynomials':['H sextic','C_star','F0','P','B','R'],'exact':True,'prefix_sha256':hashlib.sha256(raw[:48060]).hexdigest()}

report={
 'arithmetic':'exact rational; no rounding in the certificates',
 'theory_correspondence':theory_correspondence,
 'j_floor':{k:{'denominator':interval_json(v[0]),'gap':interval_json(v[1])} for k,v in j_floor().items()},
 'source_collar':cover_result,
 'arc':{str(k):{'D_outward_1e6':[outward(v[0]),outward(v[1],True)],'v_outward_1e6':[outward(v[2]),outward(v[3],True)]} for k,v in arc.items()},
 'published_scalar_barrier_reuse':barriers,
 'mirror':{'Bernstein_coefficients':[str(x) for x in fb],'derivatives':{k:interval_json(v) for k,v in mirror_derivatives.items()},'simple_margin':'1923/70000','zero_safe_identity_multipliers':['10 U j n (U+n)','-10 U j n (170b+30n+89)']},
 'branches':{'beta_abs':interval_json(beta_abs),'mu':interval_json(mu),'P_Bernstein':[str(x) for x in PB],'B_Bernstein':[str(x) for x in BBB],'R_Bernstein':[str(x) for x in RB],'ratio_low_Bernstein':[str(x) for x in bern1(QL,lo,hi)],'ratio_high_Bernstein':[str(x) for x in bern1(QH,lo,hi)],'derivative_identities':True},
 'collision_chart_quotients':quotients,
 'scope':'No effective collision h floor or slope cap; no numerical p exclusion width; no M witness or arcwise/global absence.'}
Path(__file__).with_name('certificate-results.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'complete':True,'source_collar':cover_result,'scalar_barrier':barriers,'mirror_caps':[20,25,65],'branch_derivative_identities':True,'quotients':quotients},indent=2))
