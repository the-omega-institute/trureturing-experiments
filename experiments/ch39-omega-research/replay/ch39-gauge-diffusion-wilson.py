"""Exact tensor-product Wilson tests of diffusion and unified coupling ratios."""
from fractions import Fraction as F

# Each representation is (integer 6Y, twice isospin, SU(3) Dynkin labels).
def rep(q,j,p,r):
    assert (q+3*j+2*(p+2*r))%6==0
    dim=(j+1)*(p+1)*(r+1)*(p+r+2)//2
    return dim,(F(q*q),F(j*(j+2),4),F(p*p+r*r+p*r+3*p+3*r,3))

quark=rep(1,1,1,0)
square=[rep(2,j,p,r) for j in (0,2) for p,r in ((0,1),(2,0))]
absolute=[rep(0,j,p,r) for j in (0,2) for p,r in ((0,0),(1,1))]
assert sum(dim for dim,_ in square)==sum(dim for dim,_ in absolute)==36
for i in range(3):
    assert 2*quark[1][i]==sum(dim*cas[i] for dim,cas in square+absolute)/72
unified=(F(1,60),F(1),F(1))
def rate(R):return sum(x*y for x,y in zip(R[1],unified))
color,weak,electron=rep(0,0,1,1),rep(0,2,0,0),rep(6,0,0,0)
assert 2*rate(color)==3*rate(weak)==10*rate(electron)
assert rate(quark)==F(21,10)
g1,g2,g3=[F(w) for w in (".467453",".63811","1.0583")]
print("1 TeV color/weak and electron/weak rates",
      float(3*g3*g3/(2*g2*g2)),float(3*g1*g1/(10*g2*g2)))
print("unified ratios",rate(color)/rate(weak),rate(electron)/rate(weak))
# A symmetric Poisson jump at a conjugacy class where both fundamental
# traces vanish saturates the correlated color/weak jump diagnostic.
phiC,phiW=F(-1,8),F(-1,3)
assert (1-phiC)*(1-phiW)==F(3,2)
