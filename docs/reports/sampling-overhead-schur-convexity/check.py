#!/usr/bin/env python3
"""Exact checks for Xiong--Chandra--Ng--Hanzo's unrestricted SOF conjecture.

Run with Python 3 and SymPy. All mathematical inputs and comparisons are exact.
The symbolic identities cover arbitrary a,t; rational grids supplement them.
This is experimental evidence and algebra checking, not a Lean proof.
"""
from fractions import Fraction as F
from itertools import product
import sympy as S


def conv(f, g):
    assert len(f) == len(g) and len(f) & (len(f) - 1) == 0
    return [sum((f[y] * g[x ^ y] for y in range(len(f))), F(0))
            for x in range(len(f))]


def char(k, x):
    return (-1) ** ((k & x).bit_count() % 2)


def fourier(f):
    return [sum((char(k, x) * f[x] for x in range(len(f))), F(0))
            for k in range(len(f))]


def inverse(f):
    eig = fourier(f)
    assert all(e != 0 for e in eig)
    q = [sum((F(char(k, x), len(f)) / eig[k] for k in range(len(f))), F(0))
         for x in range(len(f))]
    assert conv(f, q) == [F(1)] + [F(0)] * (len(f) - 1)
    return q


def gaussian_inverse(f):
    """Independent rational elimination; does not use Fourier diagonalization."""
    size = len(f)
    rows = [[f[i ^ j] for j in range(size)] + [F(i == 0)] for i in range(size)]
    for j in range(size):
        pivot = next(i for i in range(j, size) if rows[i][j])
        rows[j], rows[pivot] = rows[pivot], rows[j]
        d = rows[j][j]
        rows[j] = [x / d for x in rows[j]]
        for i in range(size):
            if i != j:
                c = rows[i][j]
                rows[i] = [x - c * y for x, y in zip(rows[i], rows[j])]
    return [row[-1] for row in rows]


def l1(v):
    return sum(map(abs, v), F(0))


def gamma(q):
    return l1(q) ** 2 - 1


def mix(f, delta):
    return [(1 - delta) * x + delta / len(f) for x in f]


def check_majorization(p, b, difference):
    ps, bs = sorted(p, reverse=True), sorted(b, reverse=True)
    gaps = [sum(ps[:k]) - sum(bs[:k]) for k in range(1, len(p) + 1)]
    assert all(x >= 0 for x in gaps)
    assert gaps[3] == difference > 0
    assert all(x == 0 for k, x in enumerate(gaps, 1) if k != 4)


def pauli_matrices(n):
    ident = S.eye(2)
    xmat = S.Matrix([[0, 1], [1, 0]])
    zmat = S.diag(1, -1)
    ymat = S.I * xmat * zmat
    mats = []
    for label in range(4 ** n):
        local = []
        for j in range(n):
            x, z = (label >> j) & 1, (label >> (n + j)) & 1
            local.append(ymat if x and z else xmat if x else zmat if z else ident)
        mats.append(S.kronecker_product(*local) if n > 1 else local[0])
    return mats


def physical_pauli_checks():
    mats = pauli_matrices(2)
    for i, j in product(range(16), repeat=2):
        uv, w = mats[i] * mats[j], mats[i ^ j]
        phase = S.trace(w.conjugate().T * uv) / 4
        assert phase in (1, -1, S.I, -S.I)
        assert uv == phase * w
    # Section III-C: p=(.97,.01,.01,.01), PTM=(1,.96,.96,.96).
    # The source rounds 25/24, -1/96 and 33/256 in its displayed decimals.
    p = [F(97, 100)] + [F(1, 100)] * 3
    mats = pauli_matrices(1)
    ptm = S.zeros(4)
    for j in range(4):
        out = sum((S.Rational(w.numerator, w.denominator) * m * mats[j] * m.conjugate().T
                   for w, m in zip(p, mats)), S.zeros(2))
        for i in range(4):
            ptm[i, j] = S.trace(mats[i] * out) / 2
    assert ptm == S.diag(1, *([S.Rational(24, 25)] * 3))
    q = inverse(p)
    assert q == [F(33, 32)] + [F(-1, 96)] * 3
    assert l1(q) == F(17, 16) and gamma(q) == F(33, 256)
    assert 1 - S.trace(ptm) / 4 == S.Rational(3, 100)
    print('FIDELITY source Section III-C: PTM=(1,24/25,24/25,24/25), GGEP=3/100')
    print('FIDELITY inverse=(33/32,-1/96,-1/96,-1/96), L1=17/16, SOF=33/256')
    print('PHYSICAL all 256 two-qubit Pauli products agree with XOR modulo phase')


def symbolic_family():
    a, t, d, size, lp, delta_gap = S.symbols('a t d N Lp Delta', real=True)
    A, C, D = 1 - 2*a*t, 1 - 4*a, 1 - 4*a - 2*a*t
    p = [1 - 3*a, a, a, a, 0, 0, 0, 0]
    b = [1 - 3*a, a, a, a*(1-t), a*t, 0, 0, 0]
    qp = [(1-a)/C, -a/C, -a/C, -a/C, 0, 0, 0, 0]
    assert all(S.factor(x - int(i == 0)) == 0 for i, x in enumerate(conv(p, qp)))
    signs = [(1 if (x & 3) == 0 else -1) * (-1 if x & 4 else 1) for x in range(8)]
    shat, eig = fourier(signs), fourier(b)
    assert shat == [0, 0, 0, 0, -4, 4, 4, 4]
    assert all(S.factor(x-y) == 0 for x, y in zip(eig[4:], [A,C,C,D]))
    signed_pairing = sum(S.Rational(shat[k], 8) / eig[k] for k in range(8))
    Lp = (1+2*a)/C
    Delta = 8*a*a*t*(1-2*a-a*t)/(A*C*D)
    assert S.factor(signed_pairing - Lp - Delta) == 0
    assert S.factor((1 - A) - 2*a*t) == 0
    assert S.factor((1 - D) - (4*a + 2*a*t)) == 0
    # Explicit positivity: 0<a<1/6,0<t<1 imply 0<A,C,D<1 and
    # 1-2a-at > 1-3a > 1/2, hence Delta > 4a^2t > 0.
    assert S.factor((1-2*a-a*t) - (1-3*a) - a*(1-t)) == 0
    chosen_d = a*a*t*(1-6*a)
    assert S.factor(chosen_d - 2*a*a*t*(S.Rational(1,2)-3*a)) == 0
    # Uniform perturbation gap: (Delta - 6d/N)/(1-d).
    lower_b = (lp+delta_gap+d*(1-8/size))/(1-d)
    exact_p = (lp+d*(1-2/size))/(1-d)
    assert S.factor(lower_b-exact_p-(delta_gap-6*d/size)/(1-d)) == 0
    # The inverse identity uses only f*U=q*U=U, sum(f)=sum(q)=1.
    z = S.symbols('z')
    assert S.expand((1-d)*(z-d)/(1-d)+d*(1-d)/(1-d)-z) == 0
    print('SYMBOLIC eight inverse equations, Fourier sign test and Delta identity: exact')
    print('SYMBOLIC admissible d=a^2*t*(1-6a) and perturbation-gap algebra: exact')


def parameter_checks():
    cases = 0
    for n, a, t in product([2,3], [F(1,60),F(1,30),F(1,12),F(1,10),F(1,7)],
                            [F(1,5),F(1,2),F(4,5)]):
        size = 4 ** n
        p = [1-3*a,a,a,a] + [F(0)]*(size-4)
        b = [1-3*a,a,a,a*(1-t),a*t] + [F(0)]*(size-5)
        qp,qb = inverse(p),inverse(b)
        Delta = 8*a*a*t*(1-2*a-a*t)/((1-2*a*t)*(1-4*a)*(1-4*a-2*a*t))
        d = a*a*t*(1-6*a)
        assert Delta > 4*a*a*t and 0 < d < min(F(1,2)-3*a, size*Delta/6)
        pp,bb = mix(p,d),mix(b,d)
        assert sum(pp) == sum(bb) == 1 and all(x>0 for x in pp+bb)
        assert pp[0] == bb[0] and 0 < 1-pp[0] < F(1,2)
        # Q_t averages exactly labels 3 and 4; every other coordinate is fixed.
        assert bb[3] == (1-t)*pp[3]+t*pp[4]
        assert bb[4] == t*pp[3]+(1-t)*pp[4]
        assert all(pp[i] == bb[i] for i in range(size) if i not in (3,4))
        qpp,qbb = inverse(pp),inverse(bb)
        assert qpp == [(x-d/size)/(1-d) for x in qp]
        assert qbb == [(x-d/size)/(1-d) for x in qb]
        check_majorization(pp, bb, (1-d)*a*min(t,1-t))
        assert l1(qbb)-l1(qpp) >= (Delta-6*d/size)/(1-d) > 0
        assert gamma(qbb) > gamma(qpp)
        cases += 1
    print(f'FAMILY {cases} rational (n,a,t) cases: full support, inverses, strict majorization, strict SOF reversal')


def exact_interior_witness():
    p = [F(x,8000) for x in [5549,797,797,797]+[5]*12]
    r = [F(x,8000) for x in [5549,797,797,401,401]+[5]*11]
    p0 = [F(7,10),F(1,10),F(1,10),F(1,10)]+[F(0)]*12
    b0 = [F(7,10),F(1,10),F(1,10),F(1,20),F(1,20)]+[F(0)]*11
    assert mix(p0,F(1,100)) == p and mix(b0,F(1,100)) == r
    raw = inverse(b0)
    v = [377,-49,-49,-13,-29,7,7,1]
    assert raw == [F(x,252) for x in v]+[F(0)]*8
    assert all(14*v[x]+2*v[x^1]+2*v[x^2]+v[x^3]+v[x^4] == 5040*int(x==0)
               for x in range(8))
    qp,qr = inverse(p),inverse(r)
    assert qp == gaussian_inverse(p) and qr == gaussian_inverse(r)
    assert all(x>0 for x in p+r) and sum(p) == sum(r) == 1
    assert F(1)-p[0] == F(1)-r[0] == F(2451,8000) < F(1,2)
    assert l1(qp) == F(1607,792) and l1(qr) == F(3809,1782)
    assert l1(qr)-l1(qp) == F(773,7128)
    assert gamma(qr)-gamma(qp) == F(22957327,50808384)
    check_majorization(p,r,F(99,2000))
    print('WITNESS GGEP=2451/8000; L1(p^-1)=1607/792; L1(r^-1)=3809/1782')
    print('WITNESS strict prefix-4 gap=99/2000; SOF(r)-SOF(p)=22957327/50808384 > 0')
    print('WITNESS 32 convolution equations plus two independent 16x16 rational eliminations: exact')


def main():
    physical_pauli_checks()
    symbolic_family()
    parameter_checks()
    exact_interior_witness()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
