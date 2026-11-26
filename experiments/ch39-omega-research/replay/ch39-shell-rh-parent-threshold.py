"""Zero-resolvent matching of the bulk parent, with outward-rounded mass brackets."""
from fractions import Fraction as F
from flint import acb, arb, arb_series, ctx

ctx.prec = 144
ctx.cap = 4


def zero_resolvent(s):
    # Z(s) = xi'(1/2+s)/(2s xi(1/2+s)). Use analytic values at both removable points.
    if s == 0:
        z = arb_series([arb(1)/2, 1]).zeta()
        from flint import acb_series
        psi_prime = acb_series([acb(1)/4, 1]).lgamma()[2].real*2
        return (-8 + psi_prime/4 + 2*z[2]/z[0] - (z[1]/z[0])**2)/2
    if s == arb(1)/2:
        return 1 + arb.const_euler()/2 - (4*arb.pi()).log()/2
    u = s + arb(1)/2
    z = arb_series([u, 1]).zeta()
    deriv = 1/u + 1/(u-1) - arb.pi().log()/2 + (u/2).digamma()/2 + z[1]/z[0]
    return deriv/(2*s)


def mass_bracket(s):
    target = 12*zero_resolvent(s)/13
    lower, upper = arb(0), 10+s
    def gap(mass):
        return zero_resolvent((s*s+mass*mass).sqrt()) - target
    assert gap(lower)>0 and gap(upper)<0
    for i in range(62):
        mid = (lower+upper)/2
        value = gap(mid)
        if value>0:
            lower = mid
        else:
            assert value<0
            upper = mid
    assert gap(lower)>0 and gap(upper)<0
    return lower.union(upper)


def low_traces(a,b,l,r):
    return (5*a+2*b,6*a+2*b,3*(a+b)+l,3*(a+b)+r)


for f,d in [(1,0),(0,1)]:
    a,b = 2*f,2*d
    assert 32*f+20*d-4*a+2*b == 24*(f+d)
bare = low_traces(4,4,52,52)
assert bare == (28,32,76,76)
assert F(76)/(152+F(2,3)*28) == F(57,128)
q = F(12,13)
a,b = 4+28*q,4+4*q
color,family,left,right = low_traces(a,b,52,52)
assert color == left == right == F(2140,13) and family == F(2528,13)
assert max(32*color/a,20*color/b,2*color/52) == 428

if __name__ == "__main__":
    for scale in (arb(0),arb(1)/2,arb(10),arb(100),arb(1000)):
        mass = mass_bracket(scale)
        print("Q/mu =",scale.str(8),"M/mu =",mass.str(16))
    # A separate finite zero sum for the threshold at Q=0.5 mu. Use the first
    # omitted height H and the resolvent at H to bound the tail of the difference.
    s = arb(1)/2
    mass = mass_bracket(s)
    heights = [z.imag for z in acb.zeta_zeros(1,257)]
    H = heights[-1]
    partial = sum(mass**2/((h*h+s*s)*(h*h+s*s+mass**2)) for h in heights[:-1])
    tail = 2*mass**2*zero_resolvent(H)/(H*H)
    wanted = zero_resolvent(s)/13
    assert partial<wanted and wanted<partial+tail and tail<arb("0.000002")
    print("first 256 zeros, missing difference below",tail.upper().str(8))
