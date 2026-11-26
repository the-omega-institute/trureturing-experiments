"""Baryon triality with the 54-doublet crossing and its two abelian gauge kinetic coefficients."""
from fractions import Fraction as F
from flint import arb

# All fields are left-handed Weyl fields. Charges refer to ordinary SM Y and
# X, with X_Q=1. Families without a charged N^c may keep it as a neutral state.
fields=[]
for l in (0,0,-9):
    fields.extend([(6,F(1,6),1,2,3),(3,F(-2,3),-1,1,0),(3,F(1,3),-1,1,0),
                   (2,F(-1,2),l,0,1),(1,F(1),-l,0,0),(1,F(0),-l,0,0)])
assert [sum(f[3]*f[2] for f in fields),sum(f[4]*f[2] for f in fields)]==[0,0]
for yy,xx in ((2,1),(1,2),(0,1),(0,3)):
    assert sum(dim*y**yy*x**xx for dim,y,x,_,_ in fields)==0
for k in range(1,55):
    # k colorless scalar weak doublets; a Y-neutral charge-9 scalar breaks X.
    def beta(a,b):
        fermions=sum(dim*y**a*x**b for dim,y,x,_,_ in fields)
        scalars=k*2*F(1,2)**a*0**b+0**a*9**b
        return F(2,3)*fermions+F(1,3)*scalars
    assert beta(1,1)==16 and beta(0,2)==267 and beta(2,0)==F(40+k,6)
assert [i for i in range(9) if (3*i)%9==0]==[0,3,6]


def bound():
    Q=arb("4.352e13").union(arb("4.473e13"))
    inverse=arb("38.82").union(arb("38.98"))
    p=arb("2.435e18")
    # K is the inverse squared coupling matrix in charge coordinates.
    # Diagonal kinetic matching at Q gives K_YY=5/(3*g_U²).
    invU=inverse/(4*arb.pi())
    ky=5*invU/3
    t=(p/Q).log()/(8*arb.pi()**2)
    by=arb(47)/3
    rem=ky-t*by
    assert rem>0
    needed=267*t+(16*t)**2/rem
    maxg=1/needed.sqrt()
    level=needed/invU
    assert arb(".160")<maxg<arb(".162") and 12<level<13
    # A common coupling convention g_X=g_U/sqrt(kX) thus admits kX=13
    # in this integer comparison. For kX=1 the X diagonal alone excludes 2.53Q.
    failure=(8*arb.pi()**2*invU/267).exp()
    assert failure<arb("2.51")
    # Explicit extra coupling. Largest abelian norm for any field at Planck
    # is obtained by evaluating the inverse of the positive matrix K(p).
    xx=100-267*t
    yx=-16*t
    det=rem*xx-yx**2
    assert det>0
    norms=[]
    for y,x in [(y,x) for _,y,x,_,_ in fields]+[(F(1,2),0),(F(0),9)]:
        y=arb(y.numerator)/y.denominator
        norms.append((xx*y*y-2*yx*y*x+rem*x*x)/det)
    assert max(norms)<arb("2.0")
    print("Planck one-loop upper endpoint for allowed gX",(float(maxg.lower()),float(maxg.upper())),
          "required kX >", (float(level.lower()),float(level.upper())),
          "kX=1 scale ratio at most",failure.upper().str(7),
          "max charge-coupling square at gX=.10",max(norms).upper().str(7))


if __name__=="__main__":
    bound()
