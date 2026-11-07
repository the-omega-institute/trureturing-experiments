#!/usr/bin/env python3
"""Exact checks for Gehrmann–Essler, arXiv:2605.03020v2, Eq. (23).

Run with Python 3 and SymPy. No floating-point arithmetic is used.
The full-Hamiltonian check builds physical spin matrices independently of E.
Finite checks support the proposed universal argument; they are not its proof.
"""
from collections import defaultdict
from itertools import product
import sympy as s

I = s.I
ZERO = s.zeros(2)
K = s.Matrix([[I, 1], [1, -I]])


def exact(x):
    return s.simplify(s.expand(x))


def spin(n):
    """Basis k=0,...,n represents |S-k>, S=n/2; ordinary spin ladder."""
    jp = s.zeros(n + 1)
    for k in range(1, n + 1):
        jp[k - 1, k] = s.sqrt(k * (n + 1 - k))
    jm = jp.T
    return (jp + jm) / 2, (jp - jm) / (2 * I), s.diag(*[s.Rational(n, 2) - k for k in range(n + 1)])


def density(n, D):
    x, y, z = spin(n)
    p = s.eye(n + 1) - z / s.Rational(n, 2)
    kron = s.kronecker_product
    dx, dy, dz = D
    return (kron(p, p) + 4 * (dx * (kron(y, z) - kron(z, y))
            + dy * (kron(z, x) - kron(x, z)) + dz * (kron(x, y) - kron(y, x)))).applyfunc(exact)


def trace_word(word, C, B):
    m = s.eye(2)
    for k in word:
        m = m * (C if k == 0 else B if k == 1 else ZERO)
    return exact(s.trace(m))


def physical_check(n, N, D, a, label):
    dx, dy, dz = D
    delta = exact((dy - I * dx) / (dz * s.sqrt(n)))
    assert a != 0 and exact(1 - a * delta) != 0
    C = s.diag(1, exact(1 / (1 - a * delta)))
    B = K / a
    psi = {w: trace_word(w, C, B) for w in product(range(2), repeat=N)}
    psi = {w: v for w, v in psi.items() if v != 0}
    assert psi
    for w, v in psi.items():
        assert all(not (w[j] and w[(j + 1) % N]) for j in range(N))
        assert psi.get(w[1:] + w[:1], 0) == v
    h = density(n, D)
    # Columns of the full (n+1)^2 physical density, not a projected 2-level model.
    transitions = {}
    for p, q in product(range(2), repeat=2):
        transitions[p, q] = [(r, t, h[r * (n + 1) + t, p * (n + 1) + q])
                            for r, t in product(range(n + 1), repeat=2)
                            if h[r * (n + 1) + t, p * (n + 1) + q] != 0]
    out = defaultdict(lambda: s.S.Zero)
    for w, amp in psi.items():
        for j in range(N):
            k = (j + 1) % N
            for p, q, val in transitions[w[j], w[k]]:
                target = list(w)
                target[j], target[k] = p, q
                out[tuple(target)] += val * amp
    assert all(exact(v) == 0 for v in out.values())
    leakage = sum(any(k >= 2 for k in w) for w in out)
    if dx == dy == 0:
        assert psi == {(0,) * N: s.Integer(2)}
    print(f'{label}: (2S,N)=({n},{N}), full_dimension={(n+1)**N}, '
          f'nonzero_input={len(psi)}, residual_targets={len(out)}, '
          f'leakage_targets={leakage}, H_Psi=0, TI=true, blockade=true, nonzero=true')


def symbolic_local_check():
    S, r, q, d, t, c, v = s.symbols('S r q d t c v', nonzero=True)
    u = d * r * t * (c - 1) / c
    dx, dy = (v - u) / (2 * I), (u + v) / 2
    jp = s.Matrix([[0, r, 0], [0, 0, q], [0, 0, 0]])
    x, y, z = (jp + jp.T) / 2, (jp - jp.T) / (2 * I), s.diag(S, S-1, S-2)
    p = s.eye(3) - z / S
    kron = s.kronecker_product
    h = kron(p, p) + 4 * (dx * (kron(y, z) - kron(z, y))
        + dy * (kron(z, x) - kron(x, z)) + d * (kron(x, y) - kron(y, x)))
    C, B = s.diag(1, c), t * K
    A = [C, B, ZERO]
    E = [-2*S*r*v*B + 2*I*d*r*r*s.diag(c, 1),
         -2*r*u*s.diag(S+c, S*c+1), -2*S*q*u*B]
    assert B * B == ZERO
    def reduced(expr):
        numerator = s.fraction(s.cancel(expr))[0]
        return s.Poly(numerator, r).rem(s.Poly(r*r - 2*S, r)).as_expr().expand()
    entries = []
    leakage_control = []
    for a, b in product(range(3), repeat=2):
        F = sum((h[3*a+b, 3*j+k] * A[j] * A[k]
                 for j, k in product(range(3), repeat=2)), ZERO)
        R = F - E[a] * A[b] + A[a] * E[b]
        entries.extend(reduced(x) for x in R)
        if a == 0 and b == 2:
            leakage_control.extend(reduced(x) for x in F)
    assert len(entries) == 36 and all(x == 0 for x in entries)
    assert any(x != 0 for x in leakage_control)  # dropping E_2 fails
    print('SYMBOLIC_LOCAL: 36/36 scalar residuals zero modulo r^2=2S; q symbolic; omitted E_2 rejected')


def gap_formula(word, delta, c):
    N = len(word)
    hits = [j for j, v in enumerate(word) if v]
    if not hits:
        return 1 + c**N
    if any(word[j] and word[(j + 1) % N] for j in range(N)):
        return s.S.Zero
    gaps = [(hits[(j + 1) % len(hits)] - hits[j] - 1) % N for j in range(len(hits))]
    return s.expand((-I * delta)**len(hits) * c**len(hits)
                    * s.prod(sum(c**k for k in range(g)) for g in gaps))


def polynomial_checks():
    c, t = s.symbols('c t')
    delta = s.Integer(2) + I
    for N in (2, 4, 6, 8):
        polynomials = []
        for w in product(range(2), repeat=N):
            f = gap_formula(w, delta, c)
            if N <= 6:
                direct = trace_word(w, s.diag(1, c), delta*c/(c-1)*K)
                assert s.cancel(direct - f) == 0
            assert s.cancel(c**N * f.subs(c, 1/c) - f) == 0
            polynomials.append(s.Poly(f, c))
        P = s.Matrix([[f.nth(k) for k in range(N + 1)] for f in polynomials])
        rank = P.rank()
        assert rank == N // 2 + 1
        M = N // 2
        # The source's Taylor coefficients in t=1/a, obtained from the exact rational map.
        conversion = s.Matrix([[s.expand(s.series((t/(t-delta))**j, t, 0, M+1).removeO()).coeff(t,k)
                                for j in range(N+1)] for k in range(M+1)])
        assert (P * conversion.T).rank() == M + 1
        values = list(range(2, N + 3))  # N+1 distinct admissible c, all outside {0,1}.
        evaluations = s.Matrix([[f.eval(v) for v in values] for f in polynomials])
        assert evaluations.rank() == M + 1
        print(f'GAP_POLYNOMIAL: N={N}, words={2**N}, trace_crosscheck={N<=6}, '
              f'palindromic=true, coefficient_rank={rank}, first_Taylor_rank={M+1}, family_sample_rank={M+1}')


def main():
    print('ARITHMETIC: SymPy exact integers, rationals, I and radicals; no floats')
    # Fidelity: source's verified S <= 2 construction, before higher-spin tests.
    a = s.Rational(3, 5) + I*s.Rational(1, 7)
    for n in (1, 2, 3, 4):
        r = s.sqrt(n)
        physical_check(n, 4, (2*r, 3*r, s.Integer(5)), a, 'SOURCE_FIDELITY')
    symbolic_local_check()
    for n, N in [(5,2),(5,4),(6,4),(7,4),(8,4),(5,6)]:
        r = s.sqrt(n)
        physical_check(n, N, (2*r, 3*r, s.Integer(5)), a, 'EXTENSION')
    for D in [(0,s.sqrt(5),s.Integer(1)),(s.sqrt(5),0,s.Integer(-2)),(0,0,s.Integer(2))]:
        physical_check(5, 4, D, a, 'DEGENERATE_COUPLING')
    polynomial_checks()
    print('ALL_CHECKS_OK')


if __name__ == '__main__':
    main()
