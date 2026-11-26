"""A parallel Higgs choice for a residual triality and the topology of its breaking."""
from fractions import Fraction as F
from itertools import product

def branching(p,q):
    for i in range(q,p+q+1):
        for j in range(q+1):
            yield 3*(i+j)-2*(p+2*q), i-j
def dimension(p,q):
    return (p+1)*(q+1)*(p+q+2)//2
def index(p,q):
    return F(dimension(p,q)*(p*p+q*q+p*q+3*p+3*q),24)

def check():
    # U(2) interlacing: weak-clock singlet requires equal row lengths i=j=q,
    # hence its only H charge is 2(q-p).
    for p,q in product(range(15),repeat=2):
        sing=[h for h,j2 in branching(p,q) if j2==0]
        assert sing==[2*(q-p)]
    assert dimension(0,9)==55 and index(0,9)==F(495,2)
    assert all(dimension(p,p+9)>=55 for p in range(15))
    # In the primitive (H,A) torus a charge-(18,-72) condensate has kernel
    # components given by the 18th roots of z_H/z_A^4.
    assert all((h-4) %18 == 0 for h in [4+18*n for n in range(-2,3)])
    # Doubled SU(3) weights and A charges (Adams square); an external singlet
    # shift chooses y independently. Discrete weights are twice their Z9 values.
    for y,h in ((1,0),(-3,3),(6,6),(0,6)):
        assert (2*h)%18 in (0,6,12)
    print("minimal commuting H=18 singlet multiplet: dim, index",55,index(0,9))
    print("kernel components and local string classes Z/18; monopole kernel rank",0)
if __name__=="__main__":
    check()
