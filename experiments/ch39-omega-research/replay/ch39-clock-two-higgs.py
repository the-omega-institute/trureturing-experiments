"""Replace the charge-18 elementary singlet direction by two clock-singlet fundamental directions."""
from fractions import Fraction as F
from math import gcd
from flint import arb

# Higgs rows in the broken torus coordinates H and C; the unbroken charge is 4H+a.
rows=((2,-4),(2,5))
assert rows[0][0]*rows[1][1]-rows[1][0]*rows[0][1]==18
assert gcd(-4,5)==1
# Primitive order-18 generator in (theta_H-4theta_A,theta_C)/(2pi).
theta=(F(1,18),F(7,9))
assert all((h*theta[0]+c*theta[1]).denominator==1 for h,c in rows)
assert (F(7,9)*3).denominator==3
assert min(max(abs(a),abs(b)) for a in range(-9,10) for b in range(-9,10)
           if abs(b-a)==9 and gcd(a,b)==1)==5

def check():
    # Three SM families and a charge-conjugate neutrino for each.
    sums=[0]*6
    for gen in range(3):
        for dim,tc,tw,y,B,L in ((6,2,3,1,1,0),(3,1,0,-4,-1,0),(3,1,0,2,-1,0),
                               (2,0,1,-3,0,1),(1,0,0,6,0,-1),(1,0,0,0,0,-1)):
            c=B-(9*L if gen==2 else 0)
            vals=(tc*c,tw*c,dim*c,dim*c**3,dim*y*y*c,dim*y*c*c)
            sums=[s+x for s,x in zip(sums,vals)]
    assert sums==[0]*6
    pi=arb.pi()
    gy2=12*pi/(5*arb("40.6142"))
    gx2=32*gy2/3
    ga2=gy2/18
    gc2=arb(".04")
    assert 8**2*ga2+5**2*gc2+4*gx2/3<4*pi
    assert 9**2*gc2<4*pi
    a2=sum(3*8**2 for _ in rows)//3
    assert a2==128
    print("bare abelian Higgs bAA",int(a2),"maximum combined squared breaker charge",
          (64*ga2+25*gc2+4*gx2/3).upper().str(12))
    print("fermion anomalies of C=3B-9L_tau",sums)
if __name__=="__main__":
    check()
