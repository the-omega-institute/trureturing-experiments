#!/usr/bin/env python3
"""Exact algebraic checks for Luo's optimized Hermite family.

These finite checks test the source identities and the proposed leading curvature
mechanism. They do not prove an entropy asymptotic or settle Problem 10.3.
Run with Python 3 and SymPy; no floating-point arithmetic is used.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import factorial
import sympy as s


@dataclass(frozen=True)
class Q5:
    """The ordered exact field Q(sqrt(5))."""
    p: F = F(0)
    r: F = F(0)

    @staticmethod
    def cast(x):
        return x if isinstance(x, Q5) else Q5(F(x))

    def __add__(self, other):
        z = self.cast(other)
        return Q5(self.p + z.p, self.r + z.r)

    __radd__ = __add__

    def __neg__(self):
        return Q5(-self.p, -self.r)

    def __sub__(self, other):
        return self + (-self.cast(other))

    def __rsub__(self, other):
        return self.cast(other) - self

    def __mul__(self, other):
        z = self.cast(other)
        return Q5(self.p*z.p + 5*self.r*z.r, self.p*z.r + self.r*z.p)

    __rmul__ = __mul__

    def __truediv__(self, other):
        z = self.cast(other)
        den = z.p*z.p - 5*z.r*z.r
        assert den != 0
        return self * Q5(z.p/den, -z.r/den)

    def __rtruediv__(self, other):
        return self.cast(other) / self

    def __pow__(self, n):
        assert isinstance(n, int) and n >= 0
        out, base = Q5(F(1)), self
        while n:
            if n % 2:
                out = out * base
            base = base * base
            n //= 2
        return out

    def sign(self):
        p, r = self.p, self.r
        if r == 0:
            return (p > 0) - (p < 0)
        if p == 0:
            return (r > 0) - (r < 0)
        if p > 0 and r > 0:
            return 1
        if p < 0 and r < 0:
            return -1
        cmp = p*p - 5*r*r
        assert cmp != 0  # sqrt(5) is irrational
        return ((cmp > 0) - (cmp < 0)) * ((p > 0) - (p < 0))


def add_poly(p, q):
    out = dict(p)
    for k, v in q.items():
        out[k] = out.get(k, Q5()) + v
    return {k: v for k, v in out.items() if v != Q5()}


def mul_poly(p, q):
    out = {}
    for k, v in p.items():
        for l, w in q.items():
            key = tuple(x+y for x, y in zip(k, l))
            out[key] = out.get(key, Q5()) + v*w
    return out


def monomial(i, j):
    v = [0]*4
    v[i] += 1
    v[j] += 1
    return tuple(v)


def coefficient_recurrence(K):
    # d_i exp(z^T K z / 2) = sum_j K_ij z_j exp(z^T K z / 2).
    @lru_cache(None)
    def coefficient(k):
        if any(x < 0 for x in k):
            return Q5()
        if sum(k) == 0:
            return Q5(F(1))
        if sum(k) % 2:
            return Q5()
        i = next(i for i, x in enumerate(k) if x)
        ans = Q5()
        for j in range(4):
            v = list(k)
            v[i] -= 1
            v[j] -= 1
            ans += K[i][j] * coefficient(tuple(v))
        return ans / k[i]
    return coefficient


def independent_exp_coefficient(K, target):
    # Separate product of the ten exp(K_ij*z_i*z_j) factors.
    terms = {(0, 0, 0, 0): Q5(F(1))}
    for i in range(4):
        for j in range(i, 4):
            power = monomial(i, j)
            value = K[i][j] / (2 if i == j else 1)
            bound = min(target[k]//v for k, v in enumerate(power) if v)
            out = {}
            for k, c in terms.items():
                for count in range(bound+1):
                    key = tuple(k[d]+count*power[d] for d in range(4))
                    if all(key[d] <= target[d] for d in range(4)):
                        out[key] = out.get(key, Q5()) + c*value**count/factorial(count)
            terms = out
    return terms.get(target, Q5())



def wick_hermite_jet(n, d, e, c, alpha, beta):
    """Independent Gaussian-moment check, including two weight derivatives.

    After Gaussian exponential tilting, U=sqrt(2a)*W has covariance I+K_t
    in the unphased variables. Work modulo the cube of t-1/2. Wick's
    recurrence computes E[product_i H_n(U_i)] without exp(Q) coefficients.
    """
    zero = (Q5(),)*3
    one = (Q5(F(1)), Q5(), Q5())
    def add(x, y):
        return tuple(a+b for a,b in zip(x,y))
    def scale(x, a):
        return tuple(v*a for v in x)
    def mul(x, y):
        return tuple(sum((x[j]*y[k-j] for j in range(k+1)),Q5()) for k in range(3))
    V = []
    for i in range(4):
        row = []
        for j in range(4):
            if i == j:
                row.append((1-d, alpha*(1 if i%2==0 else -1), Q5()))
            elif i//2 == j//2:
                row.append((-e, Q5(), -alpha))
            elif i%2 == j%2:
                row.append((c, beta*(1 if i%2==0 else -1), Q5()))
            else:
                row.append((c, Q5(), -beta))
        V.append(row)
    @lru_cache(None)
    def moment(k):
        if sum(k) == 0:
            return one
        i = next(i for i, x in enumerate(k) if x)
        ans = zero
        for j in range(4):
            count = k[j] - (i == j)
            if count > 0:
                v = list(k)
                v[i] -= 1
                v[j] -= 1
                ans = add(ans, scale(mul(V[i][j], moment(tuple(v))), count))
        return ans
    ans = zero
    terms = [(n-2*r, F((-1)**r*factorial(n),2**r*factorial(r)*factorial(n-2*r)))
             for r in range(n//2+1)]
    for choices in product(terms,repeat=4):
        powers = tuple(v[0] for v in choices)
        coefficient = F(1)
        for _, value in choices:
            coefficient *= value
        ans = add(ans,scale(moment(powers),coefficient))
    return tuple(v/factorial(n)**4 for v in ans)


def source_fidelity():
    q = s.symbols('q')
    a = 7/(3*s.sqrt(5))
    b = (a-1)/(a+1)
    kappa = (1-a*(1-q))/(1+a*(1-q))
    tau = (1-a+a*q/2)/(1+a-a*q/2)
    polynomial = -1+3*q/2+a*a*(1-q/2)*(1-q)**2
    factor = -(7*q-4)*(7*q*q-24*q+2)/90
    assert s.simplify(polynomial-factor) == 0
    assert s.simplify(kappa**2-tau-2*a*factor/((1+a*(1-q))**2*(1+a-a*q/2))) == 0
    assert s.simplify(b-(47-21*s.sqrt(5))/2) == 0
    assert s.simplify(kappa.subs(q, s.Rational(4,7))-(3-s.sqrt(5))/2) == 0
    print('FIDELITY: Section 7 factorization, exponential-rate identity, and endpoint values OK')
    # The source's Section 4 correction matrix, for K=1, with the
    # normalizing sqrt((2j)!) restored explicitly.
    x = s.symbols('x')
    def gaussian_moment(poly):
        ans = 0
        for (n,), c in s.Poly(s.expand(poly), x).terms():
            if n % 2 == 0:
                ans += c * (s.factorial2(n-1) if n else 1)
        return s.simplify(ans)
    A = s.Matrix(2,2,lambda l,j: gaussian_moment(x**(2*l)*s.hermite_prob(2*j,x))/(2**l*s.sqrt(s.factorial(2*j))))
    assert A == s.Matrix([[1,0],[s.Rational(1,2),1/s.sqrt(2)]])
    assert A.det() != 0
    mass, variance = s.symbols('mass variance')
    correction = A.inv()*s.Matrix([mass,variance])
    assert (A*correction - s.Matrix([mass,variance])).applyfunc(s.simplify) == s.zeros(2,1)
    print('FIDELITY: exact mass/variance correction matrix and cancellation OK')


def curvature_checks():
    a = Q5(F(0), F(7,15))
    b = (a-1)/(a+1)
    cases = [F(4,7), F(29,50), F(7,12), F(41,70), F(3,5), F(3,4)]
    comparisons = 0
    for q in cases:
        lp = (a*(1+q)-1)/(a*(1+q)+1)
        lm = (a*(1-q)-1)/(a*(1-q)+1)
        alpha, beta, kappa = (lp+lm)/2-b, (lp-lm)/2, -lm
        d,e,c = -b-alpha/2, -alpha/2, beta/2
        assert min(v.sign() for v in [d,e,c,kappa]) == 1
        assert d+e+2*c == kappa
        assert (kappa*kappa-b).sign() == 1
        tau = (1-a+a*q/2)/(1+a-a*q/2)
        assert (kappa*kappa-tau).sign() == (0 if q == F(4,7) else 1)
        K = [[(d if i==j else e) if i//2==j//2 else c for j in range(4)] for i in range(4)]
        coeff = coefficient_recurrence(K)
        # Exact midpoint derivatives of Q_t, using (sqrt(t*(1-t)))''=-2.
        first = {monomial(i,i): (-alpha/2 if i%2==0 else alpha/2) for i in range(4)}
        first = add_poly(first,{monomial(0,2):beta,monomial(1,3):-beta})
        second = {monomial(0,1):2*alpha, monomial(2,3):2*alpha,
                  monomial(0,3):-2*beta, monomial(1,2):-2*beta}
        kernel = add_poly(second,mul_poly(first,first))
        assert sum(first.values(),Q5()) == Q5()
        assert sum(second.values(),Q5()) == -4*(kappa+b)
        for n in [2,4,6,8]:
            target = (n,)*4
            C = coeff(target)
            C1 = sum((v*coeff(tuple(n-x for x in k)) for k,v in first.items()),Q5())
            C2 = sum((v*coeff(tuple(n-x for x in k)) for k,v in kernel.items()),Q5())
            assert C.sign() == 1
            assert C1 == Q5()
            assert C2.sign() == -1
            if n <= 4:
                assert C == independent_exp_coefficient(K,target)
                wick = wick_hermite_jet(n,d,e,c,alpha,beta)
                assert wick[0] == C
                assert wick[1] == C1
                assert 2*wick[2] == C2
                comparisons += 1
        # The projected exponent support includes 0 and each coordinate vector,
        # and its probability-weighted mean is exactly (1/2,1/2,1/2).
        support = []
        for i in range(4):
            for j in range(i,4):
                exponent = monomial(i,j)[:3]
                weight = K[i][j]/(2 if i==j else 1)/(2*kappa)
                assert weight.sign() == 1
                support.append((exponent,weight))
        assert sum((w for _,w in support),Q5()) == Q5(F(1))
        for j in range(3):
            assert sum((w*v[j] for v,w in support),Q5()) == Q5(F(1,2))
        assert {(0,0,0),(1,0,0),(0,1,0),(0,0,1)} <= {v for v,_ in support}
        print('CANDIDATE: q='+str(q)+'; n=2,4,6,8: C>0, Cprime=0, Csecond<0; tilt support/mean OK')
    print('INDEPENDENT COEFFICIENT CROSS-CHECKS:',comparisons)
    print('INDEPENDENT GAUSSIAN WICK CHECKS (value and two derivatives):',comparisons)
    # Separate symbolic chain rule for the norm identity J''=2<M,M''>
    # at M'=0, and the envelope-contact upper bound for its exponential rate.
    t,k,b = s.symbols('t k b',positive=True)
    kt = k-(k+b)*(s.Rational(1,2)-s.sqrt(t*(1-t)))
    assert s.simplify(s.diff(s.log(2*kt),t,2).subs(t,s.Rational(1,2))+2*(k+b)/k) == 0
    print('CANDIDATE: exact log total-mass curvature = -2*(kappa+b)/kappa OK')
    print('LIMITATION: no finite check proves the C2 local-limit remainder, entropy differentiation, or eventual entropy sign')


if __name__ == '__main__':
    source_fidelity()
    curvature_checks()
    print('ALL_CHECKS_OK')
