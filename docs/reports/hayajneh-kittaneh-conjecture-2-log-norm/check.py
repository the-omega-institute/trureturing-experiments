#!/usr/bin/env python3
"""Exact checks for MIA 21 (2018), Conjecture 2 and its Theorem 3.

SymPy checks polynomial identities; Fraction interval arithmetic encloses
exponentials using Taylor series with rigorous remainder estimates.
No floating-point number is used as evidence. This is not a Lean proof.
"""
from dataclasses import dataclass
from fractions import Fraction as Q
from math import factorial
import sympy as s


@dataclass(frozen=True)
class I:
    lo: Q
    hi: Q

    def __post_init__(self):
        assert self.lo <= self.hi

    @staticmethod
    def lift(x):
        return x if isinstance(x, I) else I(Q(x), Q(x))

    def __add__(self, other):
        z = self.lift(other)
        return I(self.lo + z.lo, self.hi + z.hi)

    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, other):
        return self + -self.lift(other)

    def __rsub__(self, other):
        return self.lift(other) + -self

    def __mul__(self, other):
        z = self.lift(other)
        products = [a * b for a in (self.lo, self.hi) for b in (z.lo, z.hi)]
        return I(min(products), max(products))

    __rmul__ = __mul__

    def square(self):
        if self.lo >= 0:
            return I(self.lo**2, self.hi**2)
        if self.hi <= 0:
            return I(self.hi**2, self.lo**2)
        return I(Q(0), max(self.lo**2, self.hi**2))


def exp_q(q):
    """Enclose exp(q), -1 <= q <= 1, by exact Taylor bounds.

    Negative q: alternating series, even/odd partial sums at 20/21.
    Positive q: positive tail bounded geometrically, since successive
    ratios from the term of degree 21 onward are at most q/22 < 1.
    """
    q = Q(q)
    assert -1 <= q <= 1
    partial = sum((q**k / factorial(k) for k in range(21)), Q(0))
    next_term = q**21 / factorial(21)
    if q < 0:
        return I(partial + next_term, partial)
    return I(partial, partial + next_term / (1 - q / 22))


def exp_i(q):
    q = I.lift(q)
    return I(exp_q(q.lo).lo, exp_q(q.hi).hi)


def frob2(T):
    return sum((z**2 for z in T), s.Integer(0))


def trace2_negative_det(T):
    return (T[0, 0] - T[1, 1])**2 + (T[0, 1] + T[1, 0])**2


def symbolic_checks():
    u, v, b = s.symbols('u v b', real=True)
    c = v - u
    B = s.Matrix([[1, b], [b, 1]])
    At, A1t, logA = s.diag(u, 1), s.diag(v, 1), s.diag(-1, 0)
    L = At * B + B * A1t * logA
    R = At * logA * B + B * A1t
    assert (L - s.Matrix([[-c, u*b], [(1-v)*b, 1]])).applyfunc(s.expand) == s.zeros(2)
    assert (R - s.Matrix([[c, (1-u)*b], [v*b, 1]])).applyfunc(s.expand) == s.zeros(2)
    assert s.expand(L.det() + c + u*(1-v)*b**2) == 0
    assert s.expand(R.det() - c + v*(1-u)*b**2) == 0
    assert s.expand(v*(1-u) - c - u*(1-v)) == 0
    assert s.expand(frob2(R) - frob2(L) - 2*c*b**2) == 0
    # Fidelity: source Theorem 3 holds on the whole family because v >= u.
    # The source also states equality at t=1/2; then u=v and R=L^T.
    assert R.subs(v, u) == L.subs(v, u).T
    assert s.expand(trace2_negative_det(L) - trace2_negative_det(R)
                    - 4*c*(1-b**2)) == 0
    a, d, f, g = s.symbols('a d f g', real=True)
    T = s.Matrix([[a, d], [f, g]])
    assert s.expand(frob2(T) - 2*T.det() - trace2_negative_det(T)) == 0
    print('FIDELITY: Theorem 3 on the family: F(R)^2-F(L)^2 = 2*(v-u)*b^2 >= 0')
    print('FIDELITY: t=1/2 gives R=L^T, hence equality in every unitarily invariant norm')
    print('FAMILY: determinant formulas and trace-norm squared gap = 4*(v-u)*(1-b^2) verified')
    print('FAMILY: v*(1-u)-(v-u) = u*(1-v) > 0, so the b^2 interval is nonempty')


def det(T):
    return T[0][0]*T[1][1] - T[0][1]*T[1][0]


def trace2(T):
    return (T[0][0]-T[1][1]).square() + (T[0][1]+T[1][0]).square()


def certified_interval_checks():
    b, t = Q(7, 8), Q(3, 4)
    x = exp_q(-Q(1, 4))
    assert Q(389, 500) < x.lo <= x.hi < Q(39, 50)
    assert 1-b > 0  # eigenvalues 1-b and 1+b of B are positive
    # Uniform enclosure includes eps=0 and every 0<eps<=1/10000.
    eps = I(Q(0), Q(1, 10000))
    alpha, beta = -1+eps, -eps
    u0, u1 = exp_i(t*alpha), exp_i(t*beta)
    v0, v1 = exp_i((1-t)*alpha), exp_i((1-t)*beta)
    L = [[u0+v0*alpha, b*(u0+v1*beta)],
         [b*(u1+v0*alpha), u1+v1*beta]]
    R = [[u0*alpha+v0, b*(u0*alpha+v1)],
         [b*(u1*beta+v0), u1*beta+v1]]
    assert det(L).hi < -Q(3, 8)
    assert det(R).hi < -Q(1, 200)
    gap = trace2(L)-trace2(R)
    assert gap.lo > Q(1, 4)
    # Positive eigenvalues exp(alpha), exp(beta) lie in [exp(-1),1].
    assert -1 <= alpha.lo <= alpha.hi <= 0
    assert -1 <= beta.lo <= beta.hi <= 0
    print('EXACT EXP BOUND: 389/500 < exp(-1/4) < 39/50')
    print('WITNESS: t=3/4, b=7/8, A_eps=diag(exp(-1+eps),exp(-eps))')
    print('UNIFORM 0<=eps<=1/10000: det(L)<-3/8; det(R)<-1/200; trace-norm squared gap>1/4')
    print('DOMAIN: B has eigenvalues 1/8,15/8; for 0<eps<=1/10000 both A eigenvalues lie strictly in (exp(-1),1)')
    print('BOUNDARY eps=0: continuous (t,b) family and concrete source-admissible counterexample certified')


def main():
    symbolic_checks()
    certified_interval_checks()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
