"""Guinand-Weil evaluation of six windowed Wilson readouts, using compact Fourier support."""
from flint import acb, arb, ctx
from fractions import Fraction as F

ctx.prec = 112
pi = arb.pi()
a = arb(1)/128
B = 4*a
f0 = 1/(3*a)
def density(u, inner):
    if inner:
        return f0-u**2/(8*a**3)+u**3/(32*a**4)
    return (B-u)**3/(96*a**4)
def integrate(pieces):
    total = acb(0)
    for inner,lo,hi in ((True,arb(0),2*a),(False,2*a,B)):
        total += acb.integral(lambda z,analytic: pieces(z,inner),lo,hi,
                              abs_tol=arb("1e-23"),rel_tol=arb("1e-23"))
    assert total.imag.contains(0) and total.real.is_finite()
    return total.real
def denominator(T):
    def integrand(u,inner):
        f = density(u,inner)
        diff = u/(8*a**3)-u**2/(32*a**4) if inner else (f0-f)/u
        # (f(0)-f(u)cos(Tu))/u, and exp(-u/2)*2u/(1-exp(-2u)).
        diff += f*T**2*u*(T*u/2).sinc()**2/2
        return 4*f*(u/2).cosh()*(T*u).cos()+diff*(u/2).exp()/(1j*u).sinc()
    w = (-B/2).exp()
    return (integrate(integrand) + f0*((arb(1)/4).digamma()-pi.log())
            +2*f0*(w.atanh()+w.atan()))
def kernel(u):
    return (u/2).exp()+(-5*u/2).exp()/(-2*u).expm1()
def cosine(T,freq):
    def integrand(v,inner):
        return density(v,inner)*(kernel(freq+v)+kernel(freq-v))*(T*v).cos()
    main = integrate(integrand)
    prime = arb(0)
    hit = []
    assert freq>B and freq+B<arb(23).log()
    for n,p in ((2,2),(3,3),(4,2),(5,5),(7,7),(8,2),(9,3),
                (11,11),(13,13),(16,2),(17,17),(19,19)):
        gap = abs(arb(n).log()-freq)
        if gap>B:
            continue
        assert gap<B and (gap<2*a or gap>2*a)
        hit.append(n)
        prime -= arb(p).log()/arb(n).sqrt()*density(gap,gap<2*a)*(gap*T).cos()
    return main,prime,hit

def run():
    T=acb.zeta_zero(719).imag
    D=denominator(T)
    assert D>0 and D.rad()<arb("1e-13")
    print("all-zero window weight at gamma_719",D.str(17))
    for r,d0,d1,d2,c2 in ((17,6,12,14,F(81,98)),(867,244,724,804,F(8741547,12563572))):
        v,w = [(1-arb(d)/(2*r)).acos()-(1-arb(d0)/(2*r)).acos() for d in (d1,d2)]
        spec = [cosine(T,l) for l in (2*v,2*w,2*(w-v),2*(v+w))]
        assert [hits for _,_,hits in spec] == ([[],[2],[],[]] if r==17 else [[],[],[],[]])
        c=(arb(c2.numerator)/c2.denominator).sqrt()
        P,Q,M,A = [(b+p)/D for b,p,_ in spec]
        values=((1-P)/2,(1-Q)/2,(3-c**2)/4
                -(1-c**2)*(P+Q)/4-(1+c)**2*M/8-(1-c)**2*A/8)
        assert all(z.rad()<arb("1e-13") for z in values)
        print("radius squared",r,"window-average deficits",*[z.str(16) for z in values])
        if r==17:
            detune=arb(2).log()-2*w
            print("log 2 detuning and period",detune.str(15),(2*pi/detune).str(15))
            print("prime contribution to normalized second deficit",(-spec[1][1]/(2*D)).str(15))
        print("archimedean and pole terms of four normalized harmonic sums",
              *[(v/D).str(9) for v,_,_ in spec])
if __name__=="__main__":
    run()
