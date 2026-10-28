#!/usr/bin/env python3
"""Exact finite checks for almost-dimension discrimination tightness.

Python 3 + SymPy. Cyclotomic polynomial remainders certify Fourier identities;
all parameters, roots and radicals are exact. No numerical optimizer is used.
"""
import sys
import sympy as s


def equal(left, right, label):
    difference = s.simplify(s.expand(left - right))
    if (any(x != 0 for x in difference) if isinstance(difference, s.MatrixBase) else difference != 0):
        raise AssertionError((label, left, right, difference))


def source_fidelity():
    # Quantum 9, 1637, Eqs. (6)-(7): eigenvalues of the equiangular
    # Gram matrix yield the vacuum-component expression via the PGM.
    count = 0
    for n in range(2, 10):
        for q in [s.Rational(0), s.Rational(1, 4), s.Rational(1, 2),
                  s.Rational(3, 4), s.Rational(1)]:
            w = q * s.Rational(n - 1, n)
            overlap = 1 - s.Rational(n, n - 1) * w
            G = s.ones(n) * overlap + s.eye(n) * (1 - overlap)
            equal(G.trace(), n, 'source Gram normalization')
            parallel = 1 + (n - 1) * overlap
            transverse = 1 - overlap
            assert parallel >= 0 and transverse >= 0
            equal(G * s.ones(n, 1), parallel * s.ones(n, 1), 'source Gram eigenvector')
            for j in range(1, n):
                v = s.eye(n)[:, 0] - s.eye(n)[:, j]
                equal(G*v, transverse*v, 'source transverse Gram eigenvector')
            pgm = (s.sqrt(parallel) + (n - 1) * s.sqrt(transverse))**2 / n**2
            eq7 = (s.sqrt(w * (n - 1)) + s.sqrt(1 - w))**2 / n
            equal(pgm, eq7, 'published vacuum bound (7)')
            a = s.Rational(1, n)
            eq8 = a + (1 - 2*a)*w + 2*s.sqrt(a*(1-a)*w*(1-w))
            equal(eq7, eq8, 'published d=1 reduction (8) to (7)')
            count += 1
    print(f'SOURCE_FIDELITY: {count} exact vacuum/PGM cases, n=2..9; Eqs. (6)-(8) agree')


def family_checks():
    z = s.Symbol('z')
    cases = states = fourier_entries = 0
    for n in range(2, 10):
        modulus = s.Poly(s.cyclotomic_poly(n, z), z, domain=s.QQ)

        def root_sum(exponents):
            polynomial = s.Poly(sum(z**(e % n) for e in exponents), z, domain=s.QQ)
            return polynomial.rem(modulus).as_expr()

        # Both U* U and U U* are checked independently, including every entry.
        for x in range(n):
            for y in range(n):
                equal(root_sum(j*(y-x) for j in range(n))/n,
                      s.Integer(x == y), 'Fourier column Gram')
                equal(root_sum(k*(x-y) for k in range(n))/n,
                      s.Integer(x == y), 'sum of rank-one measurement effects')
                fourier_entries += 2
        # Every power used as a coordinate has exact unit modulus, with
        # conjugation z -> z^(-1); z is the primitive exp(2 pi i / n).
        for j in range(n):
            equal(root_sum([j + (-j)]), 1, 'root norm square')
        for d in range(1, n):
            a = s.Rational(d, n)
            threshold = 1 - a
            epsilons = sorted(set([q * threshold for q in
                                   [s.Rational(0), s.Rational(1, 4),
                                    s.Rational(1, 2), s.Rational(3, 4),
                                    s.Rational(1)]] +
                                  [(1 + threshold)/2, s.Rational(1)]))
            for eps in epsilons:
                saturated = eps <= threshold
                # psi_x[j] = coeff[j] z^(jx). N_x = f_x f_x*,
                # f_x[j] = z^(jx)/sqrt(n). Outer products are PSD by construction.
                coeff = ([s.sqrt((1-eps)/d)] * d +
                         [s.sqrt(eps/(n-d))] * (n-d)) if saturated else [1/s.sqrt(n)]*n
                for x in range(n):
                    norm = sum(c*c*root_sum([j*x-j*x]) for j,c in enumerate(coeff))
                    projected = sum(coeff[j]**2*root_sum([j*x-j*x]) for j in range(d))
                    equal(norm, 1, 'pure state normalization')
                    assert s.simplify(projected - (1-eps)) >= 0
                    if saturated:
                        equal(projected, 1-eps, 'exact per-state constraint')
                    amplitude = sum(coeff[j]*root_sum([j*x-j*x])/s.sqrt(n) for j in range(n))
                    if saturated:
                        F = (s.sqrt(d*(1-eps)) + s.sqrt((n-d)*eps))**2/n
                        source = a+(1-2*a)*eps+2*s.sqrt(a*(1-a)*eps*(1-eps))
                        equal(F, source, 'source Eq. (221) identity')
                        equal(amplitude**2, F, 'attaining POVM success')
                    else:
                        equal(amplitude**2, 1, 'plateau success')
                    states += 1
                success = sum((sum(coeff)/s.sqrt(n))**2 for _ in range(n))/n
                equal(success, source if saturated else 1, 'uniform averaging')
                cases += 1
    # Exact nondegenerate example stated in the draft, independent direct matrices.
    n, d, eps = 4, 2, s.Rational(1,10)
    U = s.Matrix(4, 4, lambda j,x: s.I**(j*x)/2)
    Psi = s.diag(3/s.sqrt(5),3/s.sqrt(5),1/s.sqrt(5),1/s.sqrt(5))*U
    P = s.diag(1,1,0,0)
    effects = [U[:,x]*U[:,x].conjugate().T for x in range(4)]
    equal(sum(effects,s.zeros(4)),s.eye(4),'explicit POVM completeness')
    score = 0
    for x in range(4):
        rho = Psi[:,x]*Psi[:,x].conjugate().T
        equal(rho.trace(),1,'explicit trace one')
        equal((rho*P).trace(),s.Rational(9,10),'explicit leakage')
        equal((rho*effects[x]).trace(),s.Rational(4,5),'explicit success')
        score += (rho*effects[x]).trace()/4
    equal(score,s.Rational(4,5),'explicit average')
    print(f'ATTAINMENT: {cases} parameter cases; {states} individual states; {fourier_entries} Fourier/PVM entries')
    print('RANGE: n=2..9, every 1<=d<n; leakage fractions 0,1/4,1/2,3/4,1 of threshold, and two plateau points')
    print('EXAMPLE: n=4 d=2 epsilon=1/10; each projected weight=9/10; success=4/5')
    print('SCOPE: finite exact checks support the outline; no universal theorem or Lean proof is claimed')


if __name__ == '__main__':
    print('Python',sys.version.split()[0],'SymPy',s.__version__)
    print('ARITHMETIC: exact integers, rationals, cyclotomic polynomial quotients and algebraic square roots')
    source_fidelity()
    family_checks()
    print('ALL_CHECKS_OK')
