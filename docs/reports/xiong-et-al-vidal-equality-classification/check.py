#!/usr/bin/env python3
"""Exact checks for the proposed single-scale Vidal equality classification.

No optimization oracle or floating-point numbers are used. These checks do not
establish completeness of the classification or independently compute a roof.
A matching ensemble upper bound and the source's bound (14) certify equality
for the displayed examples, conditional on that published lower bound.
"""
from itertools import combinations_with_replacement
import sympy as s

R = s.Rational

def eq(a, b):
    delta = a - b
    if isinstance(delta, s.MatrixBase):
        assert delta.applyfunc(s.simplify) == s.zeros(*delta.shape), (a, b)
    else:
        assert s.simplify(delta) == 0, (a, b)


def vec(x):
    return s.Matrix([x[i, j] for i in range(x.rows) for j in range(x.cols)])


def projector(x):
    v = vec(x)
    return v * v.conjugate().T


def fidelity(d, rho):
    v = vec(d)
    return s.simplify((v.conjugate().T * rho * v)[0])


def tail(x, k):
    gram = s.simplify(x.conjugate().T * x)
    values = []
    for value, multiplicity in gram.eigenvals().items():
        assert value.is_Rational and value >= 0, value
        values.extend([value] * multiplicity)
    values.sort(reverse=True)
    return sum(values[k:], s.Integer(0))


def boundary(a, f):
    if f <= a:
        return s.Integer(0)
    return s.simplify((s.sqrt((1-a)*f)-s.sqrt(a*(1-f)))**2)


def source_fidelity():
    # Section VI: Tr(Phi_m rho_F^iso) = F, including the stated m=4,F=0.82.
    f = s.Symbol('F', real=True)
    for m in range(2, 6):
        phi = s.eye(m) / s.sqrt(m)
        p = projector(phi)
        rho = (1-f) / (m*m-1) * (s.eye(m*m)-p) + f*p
        eq(s.trace(rho), 1)
        eq(fidelity(phi, rho), f)
    a, f0 = R(3,4), R(41,50)
    eq(boundary(a,f0), (s.sqrt(82)-3*s.sqrt(6))**2/400)
    print('SOURCE_FIDELITY: isotropic identity, m=2..5, symbolic F; m=4,F=41/50 reproduced')


def explicit_mixed_family():
    d = s.diag(2,2,1) / 3
    x = s.diag(4,2,1) / s.sqrt(21)
    y = s.diag(2,4,1) / s.sqrt(21)
    px, py = projector(x), projector(y)
    for z in [x,y]:
        eq(s.trace(projector(z)), 1)
        eq(tail(z,1),R(5,21))
        eq(fidelity(d,projector(z)),R(169,189))
    eq((vec(x).T*vec(y))[0],R(17,21))
    eq((vec(d).T*vec(x))[0],13/(3*s.sqrt(21)))
    eq(boundary(R(4,9),R(169,189)),R(5,21))
    q = s.Symbol('q',real=True)
    rho = q*px+(1-q)*py
    eq(fidelity(d,rho),R(169,189))
    eq(q*tail(x,1)+(1-q)*tail(y,1),R(5,21))
    idx = [0,4,8]
    block = rho.extract(idx,idx)
    eq(block.subs(q,R(1,2)),s.Matrix([[20,16,6],[16,20,6],[6,6,2]])/42)
    z = s.Symbol('z')
    eq(block.charpoly(z).as_expr(),z*(z*z-z+R(152,441)*q*(1-q)))
    assert rho.subs(q,R(1,2)).eigenvals() == {R(19,21):1,R(2,21):1,s.Integer(0):7}
    # The positive coefficient of q(1-q) gives rank two throughout 0<q<1.
    assert R(152,441)>0
    rot = s.Matrix([[3,1,0],[1,3,0],[0,0,1]])/s.sqrt(21)
    eq(tail(rot,1),R(5,21));eq(fidelity(d,projector(rot)),R(169,189))
    print('MIXED_FAMILY: symbolic q; F=169/189; ensemble cost=R_A(F)=5/21')
    print('MIXED_RANK: rank two for 0<q<1; half-mixture eigenvalues 19/21,2/21,0 (multiplicity 7)')
    print('ROTATED_CUTOFF: off-diagonal generator has the same fidelity and pure tail')


def finite_generator_screen():
    count = rotated = 0
    for m in range(2,5):
        for amps0 in combinations_with_replacement(range(3),m):
            amps = tuple(reversed(amps0))
            total = sum(a*a for a in amps)
            if total == 0:
                continue
            nu = [R(a*a,total) for a in amps]
            d = s.diag(*amps) / s.sqrt(total)
            for k in range(1,m):
                a = sum(nu[:k],s.Integer(0));b = 1-a
                if b == 0:
                    continue
                base = s.diag(*([1]*k+[0]*(m-k)))
                projections = [base]
                if amps[k-1] == amps[k]:
                    # Both real and complex rotations in a cutoff eigenspace.
                    for phase in [s.Integer(1),s.I]:
                        w = s.zeros(m,1);w[k-1]=R(3,5);w[k]=R(4,5)*phase
                        p = base.copy();p[k-1,k-1]=0;p += w*w.conjugate().T
                        projections.append(p);rotated += 1
                for p in projections:
                    eq(p*p,p);eq(p.conjugate().T,p);assert p.rank()==k
                    eq(s.trace(p*s.diag(*nu)),a)
                    for alpha,beta in [(2,1),(3,1),(3,2)]:
                        norm = alpha*alpha*a+beta*beta*b
                        x = s.simplify(d*(alpha*p+beta*(s.eye(m)-p))/s.sqrt(norm))
                        t = beta*beta*b/norm
                        f = (alpha*a+beta*b)**2/norm
                        assert a < f < 1
                        eq(s.trace(projector(x)),1);eq(tail(x,k),t)
                        eq(fidelity(d,projector(x)),f);eq(boundary(a,f),t)
                        # Nonuniform reference: all generators share a missing coefficient.
                        if len(set(nu))>1:
                            i,j = next((i,j) for i in range(m) for j in range(m) if nu[i]!=nu[j])
                            eq(x[i,j],0)
                        # Also check rectangular zero-padding, n=m+1.
                        padded = x.col_join(s.zeros(1,m))
                        dpad = d.col_join(s.zeros(1,m))
                        eq(tail(padded,k),t);eq(fidelity(dpad,projector(padded)),f)
                        count += 1
    print('FINITE_GENERATORS: %d exact parameter/projector/scale cases; %d rotated projections; m=2..4 and n=m,m+1' % (count,rotated))


if __name__ == '__main__':
    source_fidelity()
    explicit_mixed_family()
    finite_generator_screen()
    print('LIMIT: finite and symbolic family checks; exhaustiveness and roof attainment require the proposed proof')
    print('ALL_CHECKS_OK')
