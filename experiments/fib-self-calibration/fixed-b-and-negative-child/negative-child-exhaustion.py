from fractions import Fraction as Fr
from math import comb
from collections import defaultdict
import hashlib, json
import sympy as sp

# Original traces, late-node interpolation and zero-safe source identities.
rr,bb,uu,dd,jj,xx,hh,zz=sp.symbols('rr bb uu dd jj xx hh zz')
cc=rr*(1-rr)/bb
PP=sp.Matrix([[rr,bb],[cc,1-rr]])
JM=sp.Matrix([[0,1],[1,0]]); MM=sp.Matrix([[0,1],[1,1]])
QQ=sp.eye(2); traces=[sp.trace(QQ*PP)]
for E in (JM,MM,JM,MM,MM,MM):
    QQ=E*QQ; traces.append(sp.trace(QQ*PP))
j0=bb+cc; n0=j0+rr; T0=1+bb; U0=1+2*bb
expected=[1,j0,T0,n0,U0,U0+n0,2*U0+n0]
assert all(sp.cancel(a-b)==0 for a,b in zip(traces,expected))
assert sp.cancel(sp.trace(JM*QQ*PP)-(T0+2*j0+4*rr))==0
assert sp.cancel(sp.trace(MM*QQ*PP)-(3*U0+2*n0))==0
s=sp.Symbol('s'); hs=sp.symbols('h3:7')
Hint=sum(hs[i-3]*sp.prod((s-k)/sp.Integer(i-k)
         for k in range(3,7) if k!=i) for i in range(3,7))
rows=[(20,-45,36,-10),(10,-20,15,-4),(4,-6,4,-1),(-1,4,-6,4)]
for node,row in zip((0,1,2,7),rows):
    assert sp.expand(Hint.subs(s,node)-sum(a*b for a,b in zip(row,hs)))==0
H0=-dd*uu; H2=2*uu*(sp.Rational(147,5)-dd)/(uu+1)
B0=sp.Rational(294,5)-2*dd-H2
assert sp.cancel(dd*H2+H0*B0)==0
braces=4*H0**2*zz**2+2*H0*(dd-H0)*zz*hh-(H0+dd)*(dd+(1-2*xx)*H0)*hh**2
rank=rr**2-(uu+1)*rr/2+(uu-1)*(xx*uu-(uu-1)/2)/2
assert sp.cancel(braces.subs(zz,rr*hh/uu)-4*dd**2*hh**2*rank)==0
Phi=438272*rr**4-931328*rr**3+540800*rr**2-50728*rr-853
Lam=1+6/T0+1/U0-4/j0-4/n0
assert sp.cancel((960*j0*n0*Lam).subs(bb,sp.Rational(1,8))-Phi)==0
assert sp.discriminant(Phi,rr)==3648425367354069653778543411200
r0=Fr(1,2); b0=Fr(1,4)
jw=b0+r0*(1-r0)/b0; nw=jw+r0; Tw=1+b0; Uw=1+2*b0
lw=1+6/Tw+1/Uw-4/jw-4/nw
ew=-6*(8*b0-1)*(5*b0+2)/(5*Tw*Uw*lw); dw=24+ew
samples=[-24-ew,-ew/jw,(Fr(27,5)-ew)/Tw,-ew/nw]
pw=dw+(Uw+nw)*sum(a*b for a,b in zip((-4,15,-20,10),samples))
vw=dw+(2*Uw+nw)*sum(a*b for a,b in zip((-10,36,-45,20),samples))
xw=nw/Uw
F0w=405*xw**3+(975-18*pw+5*vw)*xw**2+(90-36*pw+5*vw)*xw-480
assert (xw,pw,vw,F0w)==(Fr(7,6),Fr(19224,515),Fr(464403,2575),Fr(5736351,4120))
print('original-world traces, interpolation, zero-drift and singular interfaces: PASS',flush=True)

# Fixed polynomial normalization; no numerical root solver is used.
x,d,p,v=sp.symbols('x d p v')
DD=x*(x+1)*(x+2)
z3=(24-d)*(x+1)*(x+2); z4=(18-d)*DD
z5=(p-d)*x*(x+2); z6=(v-d)*x*(x+1)
a0=20*z3-45*z4+36*z5-10*z6
a1=10*z3-20*z4+15*z5-4*z6
a2=4*z3-6*z4+4*z5-z6
zs=x*a1-DD*(24-d)
f=sp.Poly(5*d*DD*a2+a0*((294-10*d)*DD-5*a2),d,
          domain=sp.ZZ[x,p,v])
K=sp.Poly(4*a0**2*zs**2+2*a0*(d*DD-a0)*zs*a1
          -(a0+d*DD)*(d*DD+(1-2*x)*a0)*a1**2,d,
          domain=sp.ZZ[x,p,v])
F0=405*x**3+(975-18*p+5*v)*x**2+(90-36*p+5*v)*x-480
resultant=sp.Poly(f.resultant(K),x,p,v)
R,rem=sp.div(resultant,sp.Poly(256*DD**8*F0**2,x,p,v))
assert rem.is_zero and len(R.terms())==512 and R.degree_list()==(18,6,6)
serial=json.dumps([[list(e),int(a)] for e,a in sorted(R.terms())],
                  separators=(',',':'))
assert hashlib.sha256(serial.encode()).hexdigest()==(
    'ff6b65b00b5d251082949e2318c2a53fb3f583a86c94b3d22c3e15c9ebe038e2')
assert sp.expand(R.as_expr().subs(x,0)-10**6*(32*p-5*v-366)**2*
    (120*p*p-38*p*v-3360*p+3*v*v+522*v+23112)**2)==0
A,B,C=[f.nth(i) for i in (2,1,0)]
k0,k1,k2,k3,k4=[K.nth(i) for i in range(5)]
alpha=sp.Poly(k4*(-B**3+2*A*B*C)+k3*A*(B**2-A*C)
              -k2*A*A*B+k1*A**3,x,p,v)
zeta=sp.Poly(k4*(-B*B*C+A*C*C)+k3*A*B*C-k2*A*A*C
             +k0*A**3,x,p,v)
Q=k4*(A*A*d*d-A*B*d+B*B-A*C)+k3*A*(A*d-B)+k2*A*A
assert sp.Poly(A**3*K.as_expr()-alpha.as_expr()*d-zeta.as_expr()
               -f.as_expr()*Q,d,x,p,v).is_zero
assert len(alpha.terms())==280 and len(zeta.terms())==405
print('carrier, normalization, zero-safe drift identity: PASS',flush=True)

pc=Fr(-11); vc=Fr(-8215077071969,10**12)
delta=Fr(1,10**5); rho=Fr(1,10**4)
def ceilfr(a): return -((-a.numerator)//a.denominator)
def rnd(a): return Fr((a*10**100+Fr(1,2)).__floor__(),10**100)

# Gaussian rationals and their multiplication-operator bound.
class G:
    def __init__(self,a=0,b=0): self.a=Fr(a); self.b=Fr(b)
    def __add__(self,y):
        y=asG(y); return G(self.a+y.a,self.b+y.b)
    __radd__=__add__
    def __neg__(self): return G(-self.a,-self.b)
    def __sub__(self,y): return self+-asG(y)
    def __mul__(self,y):
        y=asG(y)
        return G(self.a*y.a-self.b*y.b,self.a*y.b+self.b*y.a)
    __rmul__=__mul__
    def inverse(self):
        a=self.a*self.a+self.b*self.b
        assert a>0
        return G(self.a/a,-self.b/a)
    def norm(self): return abs(self.a)+abs(self.b)
def asG(a): return a if isinstance(a,G) else G(a)

def base_translate(P):
    out=defaultdict(Fr)
    for (i,j,k),a in P.terms():
        for jj in range(j+1):
            for kk in range(k+1):
                out[i,jj,kk]+=int(a)*comb(j,jj)*pc**(j-jj)*comb(k,kk)*vc**(k-kk)
    return {e:a for e,a in out.items() if a}

def centre_translate(P,c):
    powers=[G(1)]
    for i in range(max(e[0] for e in P)): powers.append(powers[-1]*c)
    out={}
    for (i,j,k),a in P.items():
        for l in range(i+1):
            e=(l,j,k)
            out[e]=out.get(e,G())+a*comb(i,l)*powers[i-l]
    return out

Rbase=base_translate(R)
Abase=base_translate(alpha); Zbase=base_translate(zeta)

centre_numerators=[
    (-2566651740089873802558961611116998361064043317266727026593154, 0),
    (-2488256067754961262325253267547922037066061162427561700643890, 0),
    (-2398282646064128690056667223037546435382446711027161580913045, 0),
    (-2253614772749201210358016945669058656499587992729798853866355, 0),
    (-2252626768029364415579697482419480958426371685284207315544168, 0),
    (-2184327282123042597392334337114778818666093044664055281709777, 0),
    (-385121109334004451450104015203478788150898573648732235199962, 0),
    (327012843307482581727938913901418823375830440876425286249326, 0),
    (499330948634637916987203941699780064294575742268491893596777, 0),
    (688914854153533970918949706527175689791178031035269440088934, -99781436072170823545768407807753388134017826984639043874043),
    (688914854153533970918949706527175689791178031035269440088934, 99781436072170823545768407807753388134017826984639043874043),
    (1089365608486293207913152831232653177904114251469431075943379, 0),
    (1228632880812453051104584087024888954984266181269074315374397, 0),
    (1291164724219101032355101727757188970647736976730287977692312, 0),
    (1311181683741030341929257994409657323129246507311517404455403, 0),
    (1478357740164785651714789709792381574413887570456758534609180, 0),
    (2166512041280199954483151592637154760323376849658401764138494, 0),
    (2853192397408804809598593441466492457861449208333382059958303, 0),
]
cover_rows=[[0, 721909106, 6272996575, 1443818211], [1, 557182077, 2744986374, 1114364154], [2, 387044598, 1012121942, 774089196], [3, 165990891, 203648161630, 331981781], [4, 165209397, 202067154217, 330418793], [5, 65947751, 7174020552, 131895501], [6, 14325870, 1011764684, 28651740], [7, 15072025, 2720596422, 30144049], [8, 54594931, 1369419267, 109189861], [9, 24578373, 3938457924, None], [10, 24578373, 3938457924, None], [11, 167601053, 2119294793, 335202106], [12, 138846557, 3697194846, 277693113], [13, 5287786523, 5798853285, 10575573045], [14, 242340525, 13600452720, 484681049], [15, 7395594491, 4579822343, 14791188981], [16, 1408136604, 1715507345, 2816273208], [17, 2987593861, 1599537955, 5975187722]]
lift_rows={
    7: [[14586594191, 14586999554], [7586943746, 7595671844], [5463267752, 5467513771], [238792206797, 238812815508], [-1599970973261913641111903, -1599945587307695614140967], [-30450272111, -30424162427], [217995933763, 218052167743]],
    8: [[-3036495092, -3036431507], [544124069, 551703437], [-6137104775, -6133719207], [211213381692, 211216983401], [-8228122901268666205511242, -8228052647541463959147837], [49339735826, 49343499999], [49106062687, 49123786293]],
    11: [[39513470013, 39514629802], [-91626657, -81987651], [-88641257, -87278639], [98770035747, 98826561390], [-410691884994911143352237775977, -410685543386262384945660220184], [-13007790616, -12963733491], [1684110802767, 1684562797208]],
    12: [[36201152843, 36209078125], [757691172, 978022921], [-443634329, -435131285], [86042866710, 86429707327], [-592137278411754102768899479158, -592072187899756666311418888499], [1934924860, 2221388818], [1456588903086, 1459735383010]],
    13: [[32246783527, 32261594590], [807665996, 923094076], [157813654, 178076524], [107040468731, 107776226942], [1303956819604596619017640595028, 1304239921607934064415551651976], [-13104487955, -12574455711], [1045357858003, 1047271099185]],
    14: [[34524252484, 34593841383], [-596055426, -78395912], [-591122081, -514795708], [80183940220, 83606155556], [433383498874980159661093991820, 433827360063689078526846688806], [5406510871, 7865085933], [1314537168863, 1340341580545]],
    15: [[29571620853, 29622754453], [674743704, 1158097067], [54147348, 123344563], [100973365722, 103547335728], [-2054720351871524790332466319802341, -2052934866250343879811866504921451], [-9014939004, -7273898842], [874409293482, 877761691522]],
    16: [[25236384389, 25240768644], [-543269580, 66590], [-158465814, -152043824], [85003817129, 85245159290], [-1444619946505211843744923262684895182897, -1444493617894147092855184632939788081866], [-560171046, -426794350], [660897981930, 661610803477]],
    17: [[23340043727, 23346144545], [1072293693, 1154384911], [-273034793, -264520342], [74218684287, 74584918532], [-7820420505712346074209196724238491323742763, -7819393457382634358289515809599383839002782], [1394688884, 1567448357], [612832816789, 614196214705]],
}


sigmas={}; etas={}; qs={}
centres=[G(Fr(a,10**60),Fr(b,10**60)) for a,b in centre_numerators]
for i,c in enumerate(centres):
    tr=centre_translate(Rbase,c)
    inv=tr[1,0,0].inverse(); gamma=G(rnd(inv.a),rnd(inv.b))
    eta=sum((gamma*a).norm()*delta**(j+k)
            for (l,j,k),a in tr.items() if l==0)
    s=rho if c.b==0 else 2*rho
    q=(G(1)-gamma*tr[1,0,0]).norm()+sum(
        (gamma*a*l).norm()*s**(l-1)*delta**(j+k)
        for (l,j,k),a in tr.items() if l>=1 and (l,j,k)!=(1,0,0))
    assert eta<delta and q<Fr(3,10) and eta+q*rho<rho
    assert ceilfr(eta*10**15)==cover_rows[i][1]
    assert ceilfr(q*10**12)==cover_rows[i][2]
    if c.b==0:
        sigma=Fr(ceilfr(2*eta*10**15),10**15)
        assert sigma*10**15==cover_rows[i][3]
        assert eta/(1-q)<sigma
        sigmas[i]=sigma
    etas[i]=eta; qs[i]=q
for i,c in enumerate(centres):
    for e in centres[:i]:
        wi=rho if c.b else 0; wj=rho if e.b else 0
        assert abs(c.a-e.a)>2*rho or abs(c.b-e.b)>wi+wj
assert len(sigmas)==16
assert all(centres[i].a+rho<0 for i in range(7))
positive=[7,8,11,12,13,14,15,16,17]
assert all(centres[i].a-rho>0 for i in positive)
assert all(abs(centres[i].b)>rho for i in (9,10))
print('18 contractions, 153 separations, 16 sharpenings: PASS',flush=True)

# Integer endpoints on the outward grid 10^-80.
HH=10**80
def ceildiv(a,b): return -((-a)//b)
class I:
    def __init__(self,lo,hi=None):
        self.lo=lo; self.hi=lo if hi is None else hi
        assert self.lo<=self.hi
    @staticmethod
    def exact(a):
        a=Fr(a)*HH
        return I(a.numerator//a.denominator,ceilfr(a))
    def __add__(self,y):
        y=asI(y); return I(self.lo+y.lo,self.hi+y.hi)
    __radd__=__add__
    def __neg__(self): return I(-self.hi,-self.lo)
    def __sub__(self,y): return self+-asI(y)
    def __rsub__(self,y): return asI(y)+-self
    def __mul__(self,y):
        y=asI(y); a=[u*w for u in (self.lo,self.hi) for w in (y.lo,y.hi)]
        return I(min(a)//HH,ceildiv(max(a),HH))
    __rmul__=__mul__
    def __truediv__(self,y):
        y=asI(y); assert y.lo*y.hi>0
        a=[Fr(u*HH,w) for u in (self.lo,self.hi) for w in (y.lo,y.hi)]
        return I(min(a).__floor__(),ceilfr(max(a)))
    def __rtruediv__(self,y): return asI(y)/self
    def unit(self): return self.lo*self.hi>0
    def bounds9(self):
        return [self.lo*10**9//HH,ceildiv(self.hi*10**9,HH)]
def asI(a): return a if isinstance(a,I) else I.exact(a)
def box(a,r): return I(I.exact(a-r).lo,I.exact(a+r).hi)

def symmetric_enclosure(P,c,sigma):
    tr=centre_translate(P,c)
    a=tr[0,0,0]; assert a.b==0
    rad=sum(z.norm()*sigma**l*delta**(j+k)
            for (l,j,k),z in tr.items() if (l,j,k)!=(0,0,0))
    return box(a.a,rad)

def lift(xx,dd,pp,vv):
    h3=(24-dd)/xx; h4=18-dd
    h5=(pp-dd)/(xx+1); h6=(vv-dd)/(xx+2)
    H0=20*h3-45*h4+36*h5-10*h6
    H1=10*h3-20*h4+15*h5-4*h6
    H2=4*h3-6*h4+4*h5-h6
    H7=-h3+4*h4-6*h5+4*h6
    B0=Fr(294,5)-2*dd-H2; C0=dd*dd+B0*B0
    UU=(-dd*H0+B0*H2)/C0
    bb=(UU-1)/2; jj=UU*(24-dd)/H1; rr=xx*UU-jj
    mm=dd+(3+2*xx)*H7
    return rr,bb,mm,H1,C0

Pbox=box(pc,delta); Vbox=box(vc,delta)
lift_values={}
for i in positive:
    xx=box(centres[i].a,sigmas[i])
    aa=symmetric_enclosure(Abase,centres[i],sigmas[i])
    zz=symmetric_enclosure(Zbase,centres[i],sigmas[i])
    assert aa.unit() and xx.lo>0
    dd=-zz/aa
    rr,bb,mm,H1,C0=lift(xx,dd,Pbox,Vbox)
    assert H1.unit() and C0.lo>0
    got=[a.bounds9() for a in (dd,rr,bb,mm,aa,H1,C0)]
    assert got==lift_rows[i],(i,got)
    if i==7: assert rr.lo>HH
    elif i not in (13,15): assert bb.hi<0
    elif i==13: assert 107*HH<mm.lo and mm.hi<108*HH
    elif i==15: assert 100*HH<mm.lo and mm.hi<104*HH
    lift_values[i]=(dd,rr,bb,mm,aa,H1,C0)
print('9 positive labels, 63 lift interval pairs, units and rejections: PASS',flush=True)

# Reuse Theorem 9.1's fixed point; only sharpen its enclosure here.
seed_centre=[Fr(a,10**24) for a in (
    10504686942715102009569284,88666866113376049686925,
    5156906641684654809138004,167936210357531086294279)]
eta0=Fr(973789251791947370114393715220755890,10**60)
q0=Fr(5144072038985496620415930974177109689174603232758782855326,10**60)
seed_rad=Fr(2,10**24)
assert eta0/(1-q0)<seed_rad

def seed_response(tt,bb):
    rr=tt*bb; jj=tt+bb*(1-tt*tt); nn=jj+rr; TT=1+bb; UU=1+2*bb
    LL=bb*bb+8*bb+4; AA=(8*bb-1)*(5*bb+2)
    DD=5*(LL*jj*nn-2*TT*UU*(jj+nn))
    C0=-24*DD+3*AA*jj*nn; C1=3*AA*nn
    C2=3*(49*bb+34)*jj*nn-54*UU*(jj+nn); C3=3*AA*jj
    B5=-4*C0+15*C1-20*C2+10*C3
    B6=-10*C0+36*C1-45*C2+20*C3
    B7=-20*C0+70*C1-84*C2+35*C3
    pp=(-C0+(UU+nn)*B5)/DD; vv=(-C0+(2*UU+nn)*B6)/DD
    aa=(-C0+(TT+2*jj+4*rr)*B7)/DD
    mm=(-C0+(3*UU+2*nn)*B7)/DD
    return pp,vv,aa,mm,rr,DD,nn/UU

seed_boxes=[box(a,seed_rad) for a in seed_centre]
for w,label in ((0,15),(1,13)):
    pp,vv,aa,mm,rr,DD,xx=seed_response(*seed_boxes[2*w:2*w+2])
    assert DD.unit() and rr.lo>0 and rr.hi<HH
    assert seed_boxes[2*w+1].lo>0
    low,high=-8215077071983784,-8215077071983783
    assert low*HH<vv.lo*10**15 and vv.hi*10**15<high*HH
    assert pc-delta<Fr(-11)<pc+delta
    assert vc-delta<Fr(low,10**15)<Fr(high,10**15)<vc+delta
    assert centres[label].a-rho<Fr(xx.lo,HH)
    assert Fr(xx.hi,HH)<centres[label].a+rho
    ml,mh=((102260466,102260467) if w==0 else (107408353,107408354))
    assert ml*HH<mm.lo*10**6 and mm.hi*10**6<mh*HH
print('exact T3 fixed-point containment and labels 15,13: PASS',flush=True)
