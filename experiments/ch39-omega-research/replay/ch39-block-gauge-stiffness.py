"""Two-plaquette weak stiffness at fixed microscopic quadratic trace on the standard-model group."""
from flint import acb, arb, ctx

ctx.prec = 144
pi = arb.pi()

def blocked_stiffness(kappa, coupling):
    def weight(z):
        sin = z.sin()
        s = sin**2/2
        return sin, s, (-2*kappa*s-2*coupling*s**2).exp()
    def norm(z,analytic):
        sin,_,w = weight(z)
        return sin**2*w
    def energy(z,analytic):
        sin,s,w = weight(z)
        return (kappa+2*coupling*s)**2*sin**4*z.cos()**2*w
    denominator = acb.integral(norm,0,pi)
    numerator = acb.integral(energy,0,pi)
    assert denominator.real>0 and numerator.imag.contains(0) and denominator.imag.contains(0)
    return numerator.real/(3*denominator.real)

if __name__=='__main__':
    values = [blocked_stiffness(1,j) for j in (0,1,16,256)]
    assert all(a<b for a,b in zip(values,values[1:]))
    for j,val in zip((0,1,16,256),values):
        print('quartic coupling and blocked stiffness',j,val.str(17))
    limit = 2*arb(2).sqrt()/3*(arb(9)/4).gamma()/(arb(3)/4).gamma()
    print('asymptotic coefficient',limit.str(17))

if __name__=='__main__':
    from fractions import Fraction
    def positive_transfer(spin):
        dim=2*spin+1
        bare=Fraction(4*spin*(spin+1),3)
        prob=1/bare
        assert 1-2*prob>0
        coeff=prob**2/dim**2
        return bare*coeff/((1-prob)**2+coeff)
    assert positive_transfer(1)==Fraction(4,39)
    assert positive_transfer(2)==Fraction(4,613)
    print('positive-transfer examples with bare one',positive_transfer(1),positive_transfer(2))
