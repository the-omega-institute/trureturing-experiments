"""An integer-height proxy using shell holonomies, certified up to height 10^8 by a four-dimensional lattice."""
from fractions import Fraction as F
from flint import arb,acb,fmpz_mat,ctx

ctx.prec=256
T=acb.zeta_zero(719).imag
top=10**8
w=[]
for R,ds in ((17,(6,12,14)),(867,(244,724,804))):
    a=[(1-arb(d)/(2*R)).acos() for d in ds]
    w.extend([x-a[0] for x in a[1:]])
# Exact integer proxy; choosing it uses gamma_719.
center=F(int((T*10**60).floor().unique_fmpz()),10**60)
proxy=center.limit_denominator(10**6)
p,q=proxy.as_integer_ratio()
t0=q*T-p
energy=4*sum((t0*f).sin()**2 for f in w)
eta=(energy.sqrt()/2).asin()/arb.pi()
K=int((w[0]*(q*top+p)/arb.pi()+eta).upper().ceil().unique_fmpz())
M=int((arb(K)/eta).upper().ceil().unique_fmpz())
ratios=[f/w[0] for f in w[1:]]
A=[int((r*M+arb(1)/2).floor().unique_fmpz()) for r in ratios]
B=fmpz_mat([[1,*A],[0,M,0,0],[0,0,M,0],[0,0,0,M]])
reduced,transform=B.lll(transform=True)
assert transform*B==reduced and abs(transform.det())==1
gram=reduced*reduced.transpose()
dets=[1]+[int(fmpz_mat([[gram[i,j] for j in range(n)] for i in range(n)]).det()) for n in range(1,5)]
short2=min(F(dets[i+1],dets[i]) for i in range(4))
length2=K**2+sum((M*(1+r)*eta+arb(K)/2)**2 for r in ratios)
separation=arb(short2.numerator)/short2.denominator/length2
assert separation>1
window=(arb.pi()*eta/(w[0]*q)).upper()
assert acb.zeta_zero(718).imag<arb(p)/q-window
assert acb.zeta_zero(720).imag>arb(p)/q+window
assert abs(T-arb(p)/q)<window
assert arb(top).zeta_nzeros().unique_fmpz()==248008025
print("height proxy",p,q,"q gamma-p",t0.str(16))
print("target energy",energy.str(14),"phase tolerance",eta.str(14))
print("Gram-Schmidt minimum squared / allowed squared norm",separation.lower().str(15),
      "necessary half-width in height",window.str(14))
