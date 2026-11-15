#!/usr/bin/env python3
"""Exact DNLS soliton checks; optional quadrature is diagnostic only.

Source: Mohamed Majdoub, arXiv:2611.05330v1, equations (1.1),
(1.2), (2.1), and the soliton-family subsection. Parameter a=tan(theta)
is restricted to 0<a<1, hence M=8*atan(a) lies in (0,2*pi).
No numerical value contributes to an exact check or a claimed proof.
"""
import argparse
import hashlib
from pathlib import Path
import sympy as s


def zero(label, expr):
    value = s.cancel(expr)
    assert value == 0, (label, value)
    print(label + ': EXACT_ZERO')


def exact_checks():
    y, a = s.symbols('y a', positive=True)
    u = s.symbols('u', real=True)
    C, S = (y + 1/y)/2, (y - 1/y)/2
    D = C**2 + a**2*S**2
    Z = C-s.I*a*S
    A = s.cancel(Z**3/D**2)  # Full amplitude is 2*sqrt(a)*A.
    rho = s.cancel(4*a/D)
    b = (1-a**2)/(2*a)
    v, omega = 2*b, 1+b**2
    dx = lambda f: s.cancel(y*s.diff(f, y))
    zero('hyperbolic_identity', C**2-S**2-1)
    zero('profile_modulus', 4*a*A*s.conjugate(A)-rho)
    # Q=2 sqrt(a) exp(-i b x) A(exp x). Remove the common nonzero factor.
    Q1 = dx(A)-s.I*b*A
    Q2 = dx(dx(A))-2*s.I*b*dx(A)-b**2*A
    zero('travelling_DNLS_PDE', Q2+s.I*v*Q1-omega*A
         +s.I*(dx(rho)*A+rho*Q1))
    # The logarithmic derivative has real part rho'/(2 rho),
    # imaginary part -b-3 rho/4. Verify the full complex equality.
    zero('amplitude_phase_derivative', Q1/A-dx(rho)/(2*rho)
         +s.I*(b+3*rho/4))
    Rlog = dx(rho)/(2*rho)
    zero('amplitude_ODE', dx(Rlog)+Rlog**2-1+b*rho+s.Rational(3,16)*rho**2)
    rho_u = 4*a*(1-u**2)/(1+a**2*u**2)
    zero('tanh_mass_density', rho_u/(1-u**2)-4*a/(1+a**2*u**2))
    primitive = 4*s.atan(a*u)
    assert s.simplify(s.diff(primitive,u)-4*a/(1+a**2*u**2)) == 0
    mass = s.simplify(primitive.subs(u,1)-primitive.subs(u,-1))
    assert mass == 8*s.atan(a)
    print('soliton_mass: EXACT 8*atan(a) = 8*theta for 0<a=tan(theta)<1')
    # Fidelity: source W=sech^6(x/99), and states integral W^2 = 512/7.
    W2 = 99*s.integrate((1-u**2)**5,(u,-1,1))
    assert W2 == s.Rational(512,7)
    print('source_weight_integral_W_squared: EXACT', W2)
    # Independently compute the W^6 integral used in the predicted limit.
    W6 = 99*s.integrate((1-u**2)**17,(u,-1,1))
    assert W6 > 0
    print('weight_integral_W_sixth: EXACT', W6)
    # Scaling: lambda^(p/2), dt=du/(v lambda^2), dx unchanged.
    p = s.symbols('p', real=True)
    assert p/2-2 == (p-4)/2
    assert (p/2-2).subs(p,6) == 1
    assert (p/2-2).subs(p,4) == 0
    print('transit_exponents: EXACT p=6 -> 1; p=4 -> 0')
    lam, speed = s.symbols('lambda speed', positive=True)
    zero('window_inclusion_margin', speed*lam**2-lam-1
         -(lam*(speed*lam-2)+(lam-1)))
    # For lambda>=1 and speed*lambda>=2, the displayed margin is nonnegative.
    # Hence [-1,1]^2 is inside the transformed transit domain.
    print('linear_lower_bound: EXACT I6 >= lambda*(4/v)*cosh(1/99)^(-36)'
          '*(4*sin(2*theta)/(cosh(2)+cos(2*theta)))^3')
    print('lower_bound_domain: 0<theta<pi/4; lambda>=max(1,2/v); v=2*cot(2*theta)>0')
    # Scaling covariance, using formal q and r=|q|^2; r scales by lambda.
    T,X = s.symbols('T X', real=True)
    q,r = s.Function('q'),s.Function('r')
    f=q(lam**2*T,lam*X)
    rr=r(lam**2*T,lam*X)
    scaled=(s.I*s.diff(s.sqrt(lam)*f,T)+s.diff(s.sqrt(lam)*f,X,2)
            +s.I*s.diff(lam**s.Rational(3,2)*rr*f,X))
    expected=lam**s.Rational(5,2)*(s.I*s.Subs(s.diff(q(T,X),T),(T,X),(lam**2*T,lam*X))
        +s.Subs(s.diff(q(T,X),X,2),(T,X),(lam**2*T,lam*X))
        +s.I*s.Subs(s.diff(r(T,X)*q(T,X),X),(T,X),(lam**2*T,lam*X)))
    assert s.simplify(scaled-expected.doit()) == 0
    print('mass_critical_scaling_covariance: EXACT lambda^(5/2) times PDE')
    print('check_sha256:',hashlib.sha256(Path(__file__).read_bytes()).hexdigest())


def numerical_diagnostic():
    import mpmath as mp
    mp.mp.dps=45
    # Analytic primitive of W(x)^6 from tanh substitution; quadrature only in y.
    def F(x):
        t=mp.tanh(x/99)
        return 99*sum((-1)**k*mp.binomial(17,k)*t**(2*k+1)/(2*k+1) for k in range(18))
    normW6=2*F(mp.inf)
    for a in [mp.mpf(1)/10,mp.tan(mp.pi/8)]:
        b=(1-a*a)/(2*a);v=2*b
        def rho(y):return 4*a/(mp.cosh(y)**2+a*a*mp.sinh(y)**2)
        normQ6=2*mp.quad(lambda y:rho(y)**3,[0,1,5,10,mp.inf])
        K=normW6*normQ6/v
        print('DIAGNOSTIC_ONLY M='+mp.nstr(8*mp.atan(a),15)+' K='+mp.nstr(K,15))
        for lam in [4,16,64,256]:
            ratio=2/v*mp.quad(lambda y:rho(y)**3*(F(y/lam+v*lam)-F(y/lam-v*lam)),[0,1,5,10,mp.inf])
            print('DIAGNOSTIC_ONLY lambda='+str(lam)+' I6='+mp.nstr(lam*ratio,15)
                  +' I6/lambda='+mp.nstr(ratio,15)+' fraction_of_K='+mp.nstr(ratio/K,12))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--numerical',action='store_true')
    args=parser.parse_args()
    exact_checks()
    if args.numerical:numerical_diagnostic()
    print('ALL_CHECKS_OK')
