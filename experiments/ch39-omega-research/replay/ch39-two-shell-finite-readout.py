"""Test the two actual Behrend shell readouts at per-triangle tolerance 0.001."""
from fractions import Fraction
from pathlib import Path
import runpy
from flint import arb,acb

root=Path(__file__).resolve().parent
book=runpy.run_path(str(root/"ch39-book-shell-clock-embedding.py"))
q,pi=book["q"],book["pi"]
h=book["h"]
f=book["f"]
def angle_box(p,t,c,eps):
    a,b=p.sin()**2,t.sin()**2
    a0,a1=(a-eps).lower(),(a+eps).upper()
    b0,b1=(b-eps).lower(),(b+eps).upper()
    assert 0<a0<a1<1 and 0<b0<b1<1
    S,T=a0.union(a1),b0.union(b1)
    x=abs(S.sqrt().asin()-a.sqrt().asin()).upper()
    y=abs(T.sqrt().asin()-b.sqrt().asin()).upper()
    cross=2*c*(S*(1-S)*T*(1-T)).sqrt()
    middle=(1-S)*(1-T)+c*c*S*T
    true=(p.cos()*t.cos()+c*p.sin()*t.sin())**2
    options=[middle+cross,middle-cross]
    assert sum(not abs(z-true)>eps for z in options)==1
    return x,y

def indices(w,v,aa,bb,top):
    ratio=v/w
    eta=(bb+ratio*aa)/pi
    K=int(((w*(top+h)+aa)/pi).upper().ceil().unique_fmpz())
    base=10**35
    approx=Fraction(int((base*ratio.mid()).floor().unique_fmpz()),base).limit_denominator(K)
    a,b=approx.numerator,approx.denominator
    delta=int((b*(eta+K*abs(ratio-q(a)/b))).upper().ceil().unique_fmpz())
    inv=pow(a,-1,b)
    for r in range(-delta,delta+1):
        start=r*inv % b
        for k in range(start-((start+K)//b)*b,K+1,b):
            l=int((ratio*k+q(1)/2).floor().unique_fmpz())
            if not abs(ratio*k-l)>eta:
                yield k,l

def possibilities(top,eps):
    x,y=angle_box(h*f[2],h*f[3],book["m"]["m"]["axis"],eps)
    a,b=angle_box(h*f[0],h*f[1],book["axis"],eps)
    count=0
    both=[]
    for k,l in indices(f[2],f[3],x,y,top):
        for sign in (-1,1):
            A=sign*h+(pi*k+arb(0,x))/f[2]
            B=sign*h+(pi*l+arb(0,y))/f[3]
            if not A.overlaps(B):
                continue
            H=A.intersection(B)
            if H<0 or H>top:
                continue
            count+=1
            for sign2 in (-1,1):
                V=H-sign2*h
                pos=V*f[0]/pi
                pos2=V*f[1]/pi
                assert pos.rad()+a/pi<q(1)/2 and pos2.rad()+b/pi<q(1)/2
                # A possible integer is unique; the midpoint-nearest integer
                # is it if the interval enlarged by the allowed error contains one.
                n=int((pos.mid()+q(1)/2).floor().unique_fmpz())
                j=int((pos2.mid()+q(1)/2).floor().unique_fmpz())
                a0=sign2*h+(pi*n+arb(0,a))/f[0]
                b0=sign2*h+(pi*j+arb(0,b))/f[1]
                if H.overlaps(a0):
                    I=H.intersection(a0)
                    if I.overlaps(b0):
                        both.append(I.intersection(b0))
    return count,both

def check(top):
    one,both=possibilities(top,q(".001"))
    near=[r for r in both if r.contains(h)]
    assert len(near)==1 and near[0]>acb.zeta_zero(718).imag and near[0]<acb.zeta_zero(720).imag
    other=[r for r in both if not r.contains(h)]
    print("single-shell then double-shell candidate corridors",one,len(both),"others",
          [(float(r.mid()),float(r.rad())) for r in other[:12]],flush=True)
    pieces=0
    for H in other:
        stack=[(H,0)]
        while stack:
            seg,depth=stack.pop()
            if not acb(q(1)/2,seg).zeta().contains(0):
                pieces+=1
                continue
            assert depth<18,seg
            mid=seg.mid()
            stack.extend((a,depth+1) for a in (seg.lower().union(mid),mid.union(seg.upper())))
    print("zeta-free pieces",pieces)
if __name__=="__main__":
    check(q(10)**8)
