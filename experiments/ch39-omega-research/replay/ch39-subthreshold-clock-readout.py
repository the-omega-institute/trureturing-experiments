"""Read both clocks from the two lines below the first two-vector threshold."""
from pathlib import Path
import runpy
from flint import arb

m=runpy.run_path(str(Path(__file__).with_name("ch39-weak-clock-width.py")))
q,s=m["q"],m["s"]
a,b,c=m["z1"],m["z2"],m["z3"]
labels=q("-1.09974").union(q("-1.01335"))
def high(a,b):
    sm=a+b
    p=a*b
    return (55*sm-495-11*p)/(11*sm-55-3*p)
def r(z,s):
    return s*z/((z-1)*(z-11)-s*(2*z-11))
def line(z,rho,s):
    k=r(z,s)
    return (k+rho)**2/(q(".9")*(k*k+1)-q(".2")*k+q("1.1")*(k+1)**2*((z-1)/(z-11))**2)
def check():
    eps=q("1e-9")
    A=a+arb(0,eps)
    B=b+arb(0,eps)
    C=high(A,B)
    assert (high(a,b)-c).contains(0)
    # Centered form of the rational reconstruction, keeping correlated input
    # terms. N-C0*D vanishes at the center.
    den=11*(A+B)-55-3*A*B
    change=((55-11*c)*(A-a+B-b)+(3*c-11)*(b*(A-a)+a*(B-b)+(A-a)*(B-b)))/den
    C=c+change
    dx=abs((A-a)+(B-b)+change).upper()
    dy=abs((A-a)*B*C+a*(B-b)*C+a*b*change).upper()
    assert dx<q("1e-7") and dy<q("1e-7")
    Sum=(A+B+C-13)/2
    Prod=A*B*C/11-1-Sum
    gap=(Sum**2-4*Prod).sqrt()
    sh=(Sum+gap)/2
    sl=(Sum-gap)/2
    changes=[]
    splits=[]
    for j in range(64):
        rho=labels.lower()+(labels.upper()-labels.lower())*q(j).union(q(j+1))/64
        target=line(a,rho,s)/line(b,rho,s)
        true=line(A,rho,sh)/line(B,rho,sh)
        swapped=line(A,rho,sl)/line(B,rho,sl)
        # Keep the same geometry for both measurements. The ratio error here
        # only needs to fit in the selected tolerance, not be a sharp derivative.
        changes.append(abs(true-target).upper())
        splits.append((target-swapped).lower())
    drift=max(changes)
    split=min(splits)
    print("max change in sum and product",dx.str(12),dy.str(12))
    print("upper ratio drift, lower separation for switched clocks",drift.str(12),split.str(12))
    assert drift<q(".001")<split
if __name__=="__main__":
    check()
