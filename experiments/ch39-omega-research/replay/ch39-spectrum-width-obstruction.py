"""Portal mass floor from the scalar logarithmic width alone, at any finite multiplicity."""
from fractions import Fraction
from flint import arb,ctx
ctx.prec=160
def q(z):
    v=Fraction(z)
    return arb(v.numerator)/v.denominator
def bound(width,bins=8192):
    L=q(width)
    lh=q('.02335')
    target=2*arb.pi()**2*(lh-q('-.029827'))/3
    U=16*arb.pi()*lh
    eta=1+4*lh/(8*arb.pi()-6*lh)
    budget=U/eta-6*target/L
    assert budget>0
    def bracket(j):
        l=(-2*L*j/bins).union(-2*L*(j+1)/bins)
        A=q('78/25')*(q('3/2')*l*l+2*l+q('7/2'))
        B=q('144/5')*(3-l)
        C=24*((l-2)**2+1)
        return (l.exp()*(B**2/(4*A)-C)/(18+3*l/L)).upper()
    cost=max(q(0), max(bracket(i) for i in range(bins)))
    g2=4*arb.pi()/q('38.3316')
    w2=1+13*(-q(2)).exp()
    ferm=q('36/5')*g2**2*w2
    floor=ferm-budget*cost
    return cost,budget,floor
if __name__=='__main__':
    for width in ('3.2','3.4','3.45'):
        C,A,margin=bound(width)
        assert margin>0
        print(width,'Q/m_min ceiling',q(width).exp().str(12),'C,A,gap',
              C.str(10),A.str(10),margin.lower().str(12))
