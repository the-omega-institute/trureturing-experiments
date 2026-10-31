#!/usr/bin/env python3
"""Exact SDP certificates for arXiv:1910.00708v4, Section V.E.

SymPy algebra and rational arithmetic only. These are independent finite checks
of the proposed argument, not a Lean proof or a worldwide priority claim.
"""
from itertools import product
import sympy as S

Q = S.Rational
I2 = S.eye(2)
X = S.Matrix([[0, 1], [1, 0]])
H = S.Matrix([[1, 1], [1, -1]])


def clean(a):
    return S.simplify(S.expand_complex(a))


def eq(a, b):
    if isinstance(a, S.MatrixBase):
        assert a.shape == b.shape
        assert all(clean(t) == 0 for t in a - b), (a, b)
    else:
        assert clean(a - b) == 0, (a, b)


def kron(a, b):
    return S.kronecker_product(a, b)


def diag(a):
    return S.diag(*a.diagonal())


def tr_out(a, d, r):
    return S.Matrix(d, d, lambda i, j: sum(a[i*r+k, j*r+k] for k in range(r)))


def source_fidelity():
    count = 0
    for d in range(2, 6):
        phi = S.zeros(d*d)
        for i in range(d):
            for j in range(d):
                phi[i*d+i, j*d+j] = 1
        diagonal_correlated = S.diag(*[int(i == j) for i in range(d) for j in range(d)])
        for n in range(1, 9):
            p = Q(1, n*d)
            r = min(Q(n-1, n*d*(d-1)), p)
            alpha = p*S.eye(d*d) + r*(phi-diagonal_correlated)
            eq(diag(alpha), p*S.eye(d*d))
            eig = [p, p-r, p+(d-1)*r]
            assert all(e >= 0 and Q(1, d)-e >= 0 for e in eig)
            # Feasible interval from the three eigenspaces; optimize a linear
            # function at its rational endpoints, independently of the formula.
            lower = max(-p/Q(d-1), p-Q(1, d))
            upper = min(p, (Q(1, d)-p)/Q(d-1))
            eq(upper, r)
            for lam in (Q(0), Q(1, 4), Q(1, 2), Q(1)):
                j_dep = lam*phi + Q(1-lam, d)*S.eye(d*d)
                j_deph = lam*phi + (1-lam)*diagonal_correlated
                expected = (1-lam)/n + lam*min(Q(1), Q(d, n))
                for j in (j_dep, j_deph):
                    eq(S.trace(alpha.T*j), expected)
                    for endpoint in (lower, upper):
                        trial = p*S.eye(d*d) + endpoint*(phi-diagonal_correlated)
                        assert S.trace(trial.T*j) <= expected
                count += 1
        print(f'SOURCE_FIDELITY d={d}: depolarizing and dephasing, n=1..8, lambda=0,1/4,1/2,1 OK')
    return count


def family_certificates():
    roots = {2: -S.Integer(1), 3: (-1+S.sqrt(3)*S.I)/2,
             4: S.I, 6: (1+S.sqrt(3)*S.I)/2}
    count = 0
    for d, root in roots.items():
        v = [S.Matrix([clean(root**(j*k)) for k in range(d)]) for j in range(d)]
        ps = [(w*w.H/d).applyfunc(clean) for w in v]
        u = S.Matrix.hstack(*v)
        eq(u.H*u, d*S.eye(d))
        eq(sum(ps, S.zeros(d)), S.eye(d))
        for p in ps:
            eq(p*p, p)
            eq(diag(p), S.eye(d)/d)
        shift = S.zeros(d)
        for k in range(d):
            shift[(k+1) % d, k] = 1
        for p in ps:
            eq(shift*p*shift.H, p)
        profiles = [tuple([Q(3, 4)]+[Q(0)]*(d-1)),
                    tuple(Q((-1)**j*(j+1), d+2) for j in range(d)),
                    tuple([Q(-2, 3)]*d), tuple([Q(0)]*d)]
        if d in (2, 3):
            profiles += list(product((Q(-3, 4), Q(0), Q(3, 4)), repeat=d))
        profiles = list(dict.fromkeys(profiles))
        for c in profiles:
            rs = [(I2+t*X)/2 for t in c]
            # Input-first Choi convention: J_N = sum P_j^T tensor R_j.
            jn = sum((kron(p.T, r) for p, r in zip(ps, rs)), S.zeros(2*d))
            jt = jn.T
            eq(tr_out(jn, d, 2), S.eye(d))
            assert all((1+t)/2 > 0 and (1-t)/2 > 0 for t in c)
            basis = kron(u, H)
            eig = [s for t in c for s in ((1+t)/2, (1-t)/2)]
            eq(basis.H*jt*basis, 2*d*S.diag(*eig))
            m = max(abs(t) for t in c)
            assert all((1+m)/2-e >= 0 for e in eig)
            k = max(range(d), key=lambda j: abs(c[j]))
            sign = lambda t: S.Integer(1 if t >= 0 else -1)
            rho = ps[k]
            alpha_m = kron(rho, (I2+sign(c[k])*X)/2)
            eq(S.trace(alpha_m), 1)
            eq(diag(alpha_m), kron(diag(rho), I2)/2)
            # alpha_m and rho tensor I - alpha_m are positive tensor products
            # of the explicitly verified projector factors.
            eq(alpha_m*alpha_m, alpha_m)
            gap_m = kron(rho, I2)-alpha_m
            eq(gap_m*gap_m, gap_m)
            fm = (1+m)/2
            eq(S.trace(jt*alpha_m), fm)
            blocks = [(I2+sign(t)*X)/(2*d) for t in c]
            alpha_d = sum((kron(p, a) for p, a in zip(ps, blocks)), S.zeros(2*d))
            eq(diag(alpha_d), S.eye(2*d)/(2*d))
            eq(S.trace(alpha_d), 1)
            eq(basis.H*alpha_d*basis,
               S.diag(*[s for t in c for s in (1+sign(t), 1-sign(t))]))
            eq(basis.H*(S.eye(2*d)/d-alpha_d)*basis,
               S.diag(*[s for t in c for s in (1-sign(t), 1+sign(t))]))
            fd = (1+sum(abs(t) for t in c)/d)/2
            eq(S.trace(jt*alpha_d), fd)
            # Enumerate the vertices of the real off-diagonal box, whose
            # universal validity is certified symbolically below.
            box_max = max(Q(1, 2)+sum(t*s for t, s in zip(c, signs))/(2*d)
                          for signs in product((-1, 1), repeat=d))
            eq(box_max, fd)
            eq(kron(shift, I2)*jt*kron(shift.H, I2), jt)
            count += 1
        q = Q(3, 4)
        assert (1+q)/2 > Q(3, 4)
        assert (1+q/d)/2 < Q(3, 4)
        for n in range(3, 25):
            assert (1+q)/n < Q(3, 4)
        print(f'FAMILY d={d}: {len(profiles)} profiles, Choi normalization/spectrum, primal certificates, cyclic symmetry and exact reduced optima OK')
        print(f'SEPARATION d={d}, q=3/4: F_MISC=7/8, F_DISC={(1+q/d)/2}; epsilon=1/4; largest dimensions MISC=2, DISC=1')
    return count


def symbolic_certificates():
    a, b, z, r = S.symbols('a b z r', real=True)
    block = S.Matrix([[a, z], [z, b]])
    plus, minus = S.Matrix([1, 1]), S.Matrix([1, -1])
    # Positivity of A and rI-A makes each summand nonnegative, proving
    # -r/2 <= Re A_01 <= r/2, including complex Hermitian blocks.
    for sign, left, right in ((1, plus, minus), (-1, minus, plus)):
        eq(((left.T*block*left)[0] + (right.T*(r*I2-block)*right)[0])/2,
           r+2*sign*z)
    d, q = S.symbols('d q', positive=True)
    eq((1+q)/2-Q(3, 4), (2*q-1)/4)
    eq(Q(3, 4)-(1+q/d)/2, (d-2*q)/(4*d))
    # d>=2 and 1/2<q<1 give positive margins. All n>=3 are controlled
    # by (1+q)/n <= (1+q)/3 < 2/3, without a finite cutoff argument.
    eq(Q(2, 3)-(1+q)/3, (1-q)/3)
    print('SYMBOLIC: off-diagonal interval certificate and dimension-uniform strict margins OK')


if __name__ == '__main__':
    source_count = source_fidelity()
    family_count = family_certificates()
    symbolic_certificates()
    print(f'EXACT_CHECK_COUNTS source_parameter_triples={source_count} family_profiles={family_count}')
    print('ALL_CHECKS_OK')
