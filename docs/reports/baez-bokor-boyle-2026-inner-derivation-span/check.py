#!/usr/bin/env python3
"""Exact checks of arXiv:2607.10833v1, Proposition 8 and its span conjecture.

All algebra and ranks use SymPy over Q(i). Finite-model evidence is not a
proof of the universal conjecture. No numerical approximation is used.
"""
from functools import lru_cache
from itertools import product
import platform
import sympy as s

I, HALF = s.I, s.Rational(1, 2)


def column(x):
    return s.Matrix(list(x))


def realify(x):
    v = column(x)
    return v.applyfunc(s.re).col_join(v.applyfunc(s.im))


def hstack(cols, rows):
    return s.Matrix.hstack(*cols) if cols else s.zeros(rows, 0)


def matrix_basis(p, q):
    return [s.ImmutableMatrix(p, q, lambda r, c: int(r == a and c == b))
            for a in range(p) for b in range(q)]


class Pair:
    def __init__(self, name, bp, bm, tp, tm, kp, km):
        self.name, self.bp, self.bm = name, bp, bm
        self.tp, self.tm = lru_cache(None)(tp), lru_cache(None)(tm)
        self.kp, self.km = kp, km
        self.dp, self.dm = len(bp), len(bm)
        self.n = self.dp + self.dm
        ap = hstack([column(kp(x)) for x in bp], self.dm)
        am = hstack([column(km(y)) for y in bm], self.dp)
        self.K = s.zeros(self.n)
        self.K[:self.dp, self.dp:] = am
        self.K[self.dp:, :self.dp] = ap

    def delta(self, a, b):
        plus = hstack([column(self.tp(a, b, x)) for x in self.bp], self.dp)
        minus = hstack([-column(self.tm(b, a, y)) for y in self.bm], self.dm)
        return s.diag(plus, minus)

    def sigma(self, T):
        return self.K * T.conjugate() * self.K

    def D(self, a, b):
        return self.delta(a, b) - self.delta(self.km(b), self.kp(a))

    def fidelity(self):
        assert self.K * self.K.conjugate() == s.eye(self.n)
        outer_cases, jordan_cases = 0, 0
        for basis, other, triple, opposite in [
                (self.bp, self.bm, self.tp, self.tm),
                (self.bm, self.bp, self.tm, self.tp)]:
            for x, y, z in product(basis, other, basis):
                assert triple(x, y, z) == triple(z, y, x)
                outer_cases += 1
            # Definition 1 expanded in the two signs. The identity is complex
            # multilinear in all five arguments, so basis tuples suffice.
            for a, b, x, y, z in product(basis, other, basis, other, basis):
                lhs = triple(a, b, triple(x, y, z))
                rhs = (triple(triple(a, b, x), y, z)
                       - triple(x, opposite(b, a, y), z)
                       + triple(x, y, triple(a, b, z)))
                assert lhs == rhs
                jordan_cases += 1
        for x, y, z in product(self.bp, self.bm, self.bp):
            assert self.kp(self.tp(x, y, z)) == self.tm(
                self.kp(x), self.km(y), self.kp(z))
        for x, y, z in product(self.bm, self.bp, self.bm):
            assert self.km(self.tm(x, y, z)) == self.tp(
                self.km(x), self.kp(y), self.km(z))
        conjugation_cases = 0
        for a, b, u, v in product(self.bp, self.bm, (1, I), (1, I)):
            a, b = s.ImmutableMatrix(u*a), s.ImmutableMatrix(v*b)
            d = self.delta(a, b)
            # Proposition 8, including the minus sign on V_-.
            assert self.sigma(d) == -self.delta(self.km(b), self.kp(a))
            assert self.sigma(self.D(a, b)) == self.D(a, b)
            conjugation_cases += 1
        print(f'{self.name}: fidelity outer={outer_cases}, '
              f'Jordan={jordan_cases}, conjugation={conjugation_cases}')

    def spans(self):
        # Complex bilinearity: delta(E_i,F_j) spans inn over C.
        gs = [self.delta(a, b) for a, b in product(self.bp, self.bm)]
        ambient = self.n**2
        G = hstack([column(g) for g in gs], ambient)
        basis = [s.Matrix(self.n, self.n, g) for g in G.columnspace()]
        real_basis = basis + [I*g for g in basis]
        B = hstack([realify(g) for g in real_basis], 2*ambient)
        assert B.rank() == 2*len(basis)
        # Independently solve sigma(T)-T=0 in the full complex inner span.
        C = hstack([realify(self.sigma(g)-g) for g in real_basis], 2*ambient)
        F = B*hstack(C.nullspace(), len(real_basis))
        # D is real-bilinear. Include EVERY pair from the real bases, not samples.
        ds = [self.D(s.ImmutableMatrix(u*a), s.ImmutableMatrix(v*b))
              for a, b, u, v in product(self.bp, self.bm, (1, I), (1, I))]
        R = hstack([realify(d) for d in ds], 2*ambient)
        dim_fixed, dim_D = F.rank(), R.rank()
        joined_rank = s.Matrix.hstack(F, R).rank()
        assert joined_rank == dim_fixed == dim_D
        print(f'{self.name}: dim_C inn={len(basis)}, dim_R fixed={dim_fixed}, '
              f'dim_R span(D)={dim_D}, joined_rank={joined_rank}')
        # Symbolic coefficient absorption checks the middle slot of tm too.
        u, v = s.symbols('u v', real=True)
        c = u + I*v
        for a, b in product(self.bp, self.bm):
            diff = self.delta(s.ImmutableMatrix(c*a), b)-c*self.delta(a, b)
            assert all(s.expand(x) == 0 for x in diff)
        # Reverse-inclusion formula on an exact non-real spanning expression.
        T = sum((s.Rational(j+1, j+2)*(1+(j+2)*I)*g
                 for j, g in enumerate(gs)), s.zeros(self.n))
        A = HALF*(T+self.sigma(T))
        reconstructed = HALF*sum((self.D(
            s.ImmutableMatrix(s.Rational(j+1, j+2)*(1+(j+2)*I)*a), b)
            for j, (a, b) in enumerate(product(self.bp, self.bm))),
            s.zeros(self.n))
        assert A == reconstructed and self.sigma(A) == A
        return dim_fixed


def rectangular(p, q, sp=None, sq=None):
    sp, sq = sp or [1]*p, sq or [1]*q
    Hp, Hq = s.diag(*sp), s.diag(*sq)
    def triple(x, y, z):
        return s.ImmutableMatrix(HALF*(x*y*z+z*y*x))
    def kp(x):
        return s.ImmutableMatrix(Hq*x.conjugate().T*Hp)
    def km(y):
        return s.ImmutableMatrix(Hp*y.conjugate().T*Hq)
    tag = f'M_{p},{q}' if all(x == 1 for x in sp+sq) else f'M_{p},{q}_indefinite'
    return Pair(tag, matrix_basis(p, q), matrix_basis(q, p), triple, triple, kp, km)


def dual_numbers():
    # A=C[e]/(e^2), A^+ has Jordan product xy and pair product xyz.
    basis = [s.ImmutableMatrix([1, 0]), s.ImmutableMatrix([0, 1])]
    def mul(a, b):
        return s.ImmutableMatrix([a[0]*b[0], a[0]*b[1]+a[1]*b[0]])
    def triple(a, b, c):
        return mul(mul(a, b), c)
    return Pair('dual_numbers_nonsemisimple', basis, basis, triple, triple,
                lambda x: x.conjugate(), lambda x: x.conjugate())


def zero_pair(n):
    basis = [s.ImmutableMatrix(s.eye(n)[:, j]) for j in range(n)]
    def triple(a, b, c):
        return s.ImmutableMatrix(s.zeros(n, 1))
    return Pair(f'zero_product_dim_{n}', basis, basis, triple, triple,
                lambda x: x.conjugate(), lambda x: x.conjugate())


def main():
    print('Exact arithmetic: Q(i), symbolic real u,v; no floating point')
    print(f'Python {platform.python_version()}; SymPy {s.__version__}')
    for p, q in [(1, 1), (1, 2), (2, 1), (2, 2)]:
        model = rectangular(p, q)
        model.fidelity()
        dim = model.spans()
        if q == 1:
            # Section 4: real inner algebra for C^p is u(p).
            assert dim == p*p
            print(f'{model.name}: source u({p}) dimension={p*p}')
    model = rectangular(2, 2, [1, -1], [1, -1])
    model.fidelity()
    model.spans()
    # Hermitian Jordan algebra M_2(C)^+, involution x -> H x^dagger H.
    # Evaluate the SOURCE trace form to certify this example is nonpositive.
    x = s.ImmutableMatrix([[0, 1], [0, 0]])
    Lxx = hstack([column(model.tp(x, model.kp(x), z)) for z in model.bp], 4)
    trace_form = s.trace(Lxx)
    assert x != s.zeros(2) and trace_form == -2
    print(f'{model.name}: source trace form <E12|E12>={trace_form} < 0')
    model = dual_numbers()
    model.fidelity()
    model.spans()
    assert model.tp(model.bp[1], model.bp[1], model.bp[0]) == s.zeros(2, 1)
    for n in (0, 1, 2):
        model = zero_pair(n)
        model.fidelity()
        assert model.spans() == 0
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
