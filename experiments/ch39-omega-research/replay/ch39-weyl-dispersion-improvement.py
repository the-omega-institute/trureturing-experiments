"""Monotone quintic kinetic improvement at all three Weyl centers and the quartic obstruction."""
import sympy as sp


def verify():
    c = (sp.sqrt(5)-1)/4
    S = 1-c*c
    A,B = 5*(1-c),sp.Rational(7,2)
    alpha_max = (B-A)/(B*c)
    assert 0 < alpha_max < sp.Rational(1,20)
    assert c*c > sp.Rational(9,100)
    alpha,s = sp.symbols("alpha s",positive=True)
    d = (1-alpha+alpha*S/(2*c*c))/S**2
    G = alpha+alpha*(s-S)/(2*c*c)+d*(s-S)**2
    assert sp.simplify(G.subs(s,0)-1) == 0
    assert sp.simplify(G.subs(s,S)-alpha) == 0
    assert sp.simplify(2*c*c*sp.diff(G,s).subs(s,S)-alpha) == 0
    assert 8*sp.Rational(9,100)**2*sp.Rational(19,20) > sp.Rational(1,20)
    # A normalized cubic g(y) has G(s)=1+a*s. Its only stationary angular slope fails monotonicity.
    a = 1/(2*c*c-S)
    assert 1+a < 0
    # Compare two independent symmetric quartics on the tangent hyperplane.
    shapes = [(1,-1,0,0,0),(1,1,-1,-1,0)]
    moments = sp.Matrix([[sum(v**4 for v in l),sum(v*v for v in l)**2] for l in shapes])
    assert moments.det() != 0
    m = sp.symbols("m",positive=True)
    assert sp.simplify(c/(m-A)-1/m-(1-c)*(5-m)/(m*(m-A))) == 0
    # Expand the polar normalization, retaining all homogeneous terms through fourth order.
    t,M,p2,p4,f1,f3,cos = sp.symbols("t M p2 p4 f1 f3 cos",nonzero=True)
    norm = f1**2*t**2*p2+f1*f3*t**4*p4/3
    mass = M-cos*t**2*p2/2
    h = sp.series(2*(1-(1+norm/mass**2)**(-sp.Rational(1,2))),t,0,5).removeO()
    quartic = f1*f3*p4/(3*M**2)+(f1**2*cos/M**3-3*f1**4/(4*M**4))*p2**2
    assert sp.simplify(h.coeff(t,4)-quartic) == 0
    print("Normalized monotone quintic; fixed-Wilson quartic splitting on independent invariants.")


if __name__ == "__main__":
    verify()
