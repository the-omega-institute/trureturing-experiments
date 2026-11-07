#!/usr/bin/env python3
"""Exact checks of arXiv:2606.08227v2 and its fixed-decoder extension.

SymPy expressions, rational coefficients, and polynomial identities only.
The full tensor-product protocol supplies the probabilities; the proposed
sector formula is an independently computed comparison, not an input.
"""
import sympy as S
from functools import reduce

I, R = S.I, S.Rational
Y = S.Matrix([[0, -I], [I, 0]])
A = S.diag(I, -I)
a, x, y, h, k, c, d = S.symbols('a x y h k c d', real=True)
z = x + I*y
relations = S.groebner([k*k+h*h-1, d*d+c*c-1], k, d, h, c, a, x, y, domain=S.QQ)

def zero(expr):
    expr = S.expand(S.simplify(expr))
    assert relations.reduce(expr)[1] == 0, expr

def tensor(matrices):
    return reduce(S.kronecker_product, matrices)

def protocol(n, deviation, player):
    yn = tensor([Y]*n)
    j = (S.eye(2**n)-I*yn)/S.sqrt(2)
    local = tensor([deviation if i == player else A for i in range(n)])
    # Only the input columns are restricted; output spans the full Hilbert space.
    columns = (j.conjugate().T*local)[:, [0, 2**n-1]]
    block = S.Matrix([[a, z], [S.conjugate(z), 1-a]])
    probs = [S.expand((columns[b, :]*block*columns[b, :].conjugate().T)[0])
             for b in range(2**n)]
    def payoff(b):
        measured = [(b >> (n-1-i)) & 1 for i in range(n)]
        volunteers = [1-v for v in measured]
        count = sum(volunteers)
        return 0 if count == 0 else 2-R(1,count) if volunteers[player] else 2
    return probs, S.expand(sum(payoff(b)*p for b,p in enumerate(probs)))

# U(theta,pi/2+t), with h=cos(theta/2), k=sin(theta/2), c+i*d=exp(i*t).
u = S.Matrix([[(I*c-d)*h, k], [-k, (-I*c-d)*h]])
for n in (2,4,6):
    sign = (-1)**(n//2)
    target = S.Matrix([1/S.sqrt(2), -I*sign/S.sqrt(2)])
    block = S.Matrix([[a,z],[S.conjugate(z),1-a]])
    fidelity = S.expand((target.conjugate().T*block*target)[0])
    assert S.simplify(fidelity-(R(1,2)+sign*y)) == 0
    # Source Eq. (17): J^gamma |0>, overlap with J |0>.
    g0,g1 = S.symbols('g0 g1', real=True)
    source = S.Matrix([g0,-I*sign*g1])
    source_fidelity = S.expand((target.conjugate().T*source)[0]
                               *(source.conjugate().T*target)[0])
    assert S.expand(source_fidelity-(g0*g0+g1*g1+2*g0*g1)/2) == 0
    gamma = S.symbols('gamma', real=True)
    assert S.trigsimp(source_fidelity.subs({g0:S.cos(gamma/2),g1:S.sin(gamma/2)})
                      -(1+S.sin(gamma))/2) == 0
    B = 2-R(1,n)
    formula = h*h*B*(R(1,2)+sign*S.im((c+I*d)**2*z))+(1-h*h)*(2-fidelity)
    for player in range(n):
        probs,pay = protocol(n,u,player)
        zero(sum(probs)-1)
        zero(pay-formula)
    print(f'n={n}: Eq.(17) fidelity and full {2**n}-outcome payoff identity: OK ({n} players)')

# General optimization: after phase alignment, q ranges over [0,1].
N,F,r,q = S.symbols('N F r q', real=True)
B = 2-1/N
aligned = q*B*(R(1,2)+r)+(1-q)*(2-F)
assert S.expand(aligned-(q*aligned.subs(q,1)+(1-q)*aligned.subs(q,0))) == 0
assert S.factor(B*F-(2-F)) == (3*F*N-F-2*N)/N
print('classification algebra: endpoints B(1/2+|z|), 2-F; threshold 2n/(3n-1): OK')

# Infinite rational parameterization: t=tan(delta/2) in (0,1).
t = S.symbols('t', positive=True)
cosdelta = (1-t*t)/(1+t*t)
sindelta = 2*t/(1+t*t)
phase = cosdelta+I*sindelta
Fdelta = (1+cosdelta**2)/2
assert S.factor(1-Fdelta) == 2*t*t/(t*t+1)**2
assert S.factor(cosdelta*(1-cosdelta)/2) == -t*t*(t-1)*(t+1)/(t*t+1)**2
# These displayed factors prove positivity for 0<t<1 and convergence as t->0.
assert S.limit(Fdelta,t,0,dir='+') == 1
for n in (2,4,6):
    sign = (-1)**(n//2)
    Z = S.expand(I*sign*(1+phase**2)/4)
    assert S.simplify(S.im(Z)*sign+R(1,2)-Fdelta) == 0
    assert S.simplify(Z*S.conjugate(Z)-cosdelta**2/4) == 0
    # Eigenvalues: (1 +/- cos(delta))/2, both strictly positive.
    assert S.factor((1-cosdelta)/2) == t*t/(t*t+1)
    assert S.factor((1+cosdelta)/2) == 1/(t*t+1)
    assert S.simplify(S.im(S.conjugate(phase)*Z)*sign-cosdelta/2) == 0
    # Any nontrivial cut has the same four active product basis vectors.
    cut_pt = S.Matrix([[R(1,2),0,0,0], [0,0,Z,0],
                       [0,S.conjugate(Z),0,0], [0,0,0,R(1,2)]])
    lam = S.symbols('lambda')
    assert S.simplify(cut_pt.charpoly(lam).as_expr()
                     -(lam-R(1,2))**2*(lam**2-cosdelta**2/4)) == 0
print('family: F=(1+cos(delta)^2)/2 -> 1; gain=B*cos(delta)*(1-cos(delta))/2>0: OK')
print('family: rank 2; partial-transpose negative eigenvalue -cos(delta)/2 across every cut: OK')

# Exact two-player witness, recomputed through the complete 4x4 circuit.
z0 = R(30,169)-I*R(72,169)
rho = S.zeros(4);rho[0,0]=rho[3,3]=R(1,2);rho[0,3]=z0;rho[3,0]=S.conjugate(z0)
assert rho.eigenvals() == {R(25,26):1,R(1,26):1,S.Integer(0):2}
pt = S.Matrix(4,4,lambda i,j:rho[2*(j//2)+i%2,2*(i//2)+j%2])
assert pt.eigenvals() == {R(1,2):2,R(6,13):1,-R(6,13):1}
psi = S.Matrix([1,0,0,I])/S.sqrt(2)
assert S.simplify((psi.conjugate().T*rho*psi)[0]) == R(313,338)
assert 2*144**2 > 169**2  # F>(1+1/sqrt(2))/2, both sides positive.
v = S.diag((1+5*I)/S.sqrt(26),(1-5*I)/S.sqrt(26))
p0,g0 = protocol(2,A,0);p1,g1 = protocol(2,v,0)
sub={a:R(1,2),x:R(30,169),y:-R(72,169)}
assert [S.simplify(p.subs(sub)) for p in p0] == [R(313,338),0,0,R(25,338)]
assert [S.simplify(p.subs(sub)) for p in p1] == [R(25,26),0,0,R(1,26)]
assert S.simplify(g0.subs(sub)) == R(939,676)
assert S.simplify(g1.subs(sub)) == R(975,676)
assert S.simplify((g1-g0).subs(sub)) == R(9,169)
print('n=2 witness: eigenvalues 25/26,1/26,0,0; NPT eigenvalue -6/13; F=313/338: OK')
print('n=2 witness: payoffs 939/676 -> 975/676; strict gain 9/169: OK')
print('ALL_CHECKS_OK')
