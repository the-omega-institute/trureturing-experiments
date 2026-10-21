"""Check of the ladder construction for Conjecture 1 (iii) of arXiv:2502.02491v2 at N = 1, 2, 3 with symbolic gamma.

Operators act on polynomial wavefunctions in the position representation (2.1).  With x = p1 + i p2, y = p1 - i p2,
q+ = q1 + i q2, q- = q1 - i q2, D = q.p, C = q1 p2 - q2 p1 and Q_F(s, t) = sum_k g_k sum_j s^(k-1-j) t^j:
K+ = (x^2 + x q+ Q_F(D, i(C+2)))/2, K- = (y^2 + y q- Q_F(D, i(2-C)))/2, K = C/2, I = H/2 - (K+ + K-)/2,
I' = H/2 + (K+ + K-)/2.  The script checks (2.7), (2.9), (2.10) with (5.1) on test monomials, and reconstructs the
normal-ordered coefficients of I - p2^2 and I' - p1^2 to check the form (2.6).  It also reproduces the source's
N = 2 result: the symmetries (3.1), relation (3.2) and the ladders (3.4) satisfy (3.5), which is (5.1) at N = 2.
"""
import sympy as sp

q1, q2 = sp.symbols('q1 q2')
I = sp.I


def Q1(f): return q1 * f
def Q2(f): return q2 * f
def P1(f): return -I * sp.diff(f, q1)
def P2(f): return -I * sp.diff(f, q2)
def add(*ops): return lambda f: sum(o(f) for o in ops)


def comp(*ops):
    def g(f):
        for o in reversed(ops):
            f = o(f)
        return f
    return g


def scal(c, op): return lambda f: c * op(f)


def poly_op(coeffs, op):
    def g(f):
        res, cur = 0, f
        for k, c in enumerate(coeffs):
            if k > 0:
                cur = op(cur)
            res += c * cur
        return res
    return g


D = add(comp(Q1, P1), comp(Q2, P2))
C = add(comp(Q1, P2), scal(-1, comp(Q2, P1)))


def build(N):
    g = sp.symbols('g1:%d' % (N + 1))
    Fc = [0] + list(g)
    H = add(comp(P1, P1), comp(P2, P2), poly_op(Fc, D))
    x = add(P1, scal(I, P2)); y = add(P1, scal(-I, P2))
    qp = add(Q1, scal(I, Q2)); qm = add(Q1, scal(-I, Q2))

    def QF(sign):
        tC = (lambda f: I * (C(f) + 2 * f)) if sign > 0 else (lambda f: I * (2 * f - C(f)))

        def g_(f):
            res = 0
            for k in range(1, N + 1):
                for j in range(k):
                    h = f
                    for _ in range(j):
                        h = tC(h)
                    for _ in range(k - 1 - j):
                        h = D(h)
                    res += g[k - 1] * h
            return res
        return g_
    Kp = lambda f: (x(x(f)) + x(qp(QF(1)(f)))) / 2
    Km = lambda f: (y(y(f)) + y(qm(QF(-1)(f)))) / 2
    K = lambda f: C(f) / 2

    def Phi(shift):
        A1 = lambda f: 2 * I * (K(f) + shift * f)
        A2 = lambda f: -2 * I * (K(f) + (shift - 1) * f)
        T1 = add(H, scal(-1, poly_op(Fc, A1))); T2 = add(H, scal(-1, poly_op(Fc, A2)))
        return lambda f: T1(T2(f)) / 4
    Iop = lambda f: H(f) / 2 - (Kp(f) + Km(f)) / 2
    Ipop = lambda f: H(f) / 2 + (Kp(f) + Km(f)) / 2
    return H, Kp, Km, K, Phi, Iop, Ipop


def coefficients(L, order):
    """Normal-ordered coefficients c_(a,b)(q) of L = sum c_(a,b) d1^a d2^b, recovered from L(q1^a q2^b)."""
    c = {}
    for tot in range(order + 1):
        for a in range(tot + 1):
            b = tot - a
            mono = q1 ** a * q2 ** b
            rest = L(mono)
            for (a2, b2), cv in c.items():
                if a2 <= a and b2 <= b:
                    rest -= cv * sp.diff(mono, q1, a2, q2, b2)
            c[(a, b)] = sp.expand(rest / (sp.factorial(a) * sp.factorial(b)))
    return c


ok = True
tests = [q1 ** a * q2 ** b for a in range(4) for b in range(4) if a + b <= 4] + [q1 ** 2 * q2 ** 3]
for N in (1, 2, 3):
    H, Kp, Km, K, Phi, Iop, Ipop = build(N)
    for f in tests:
        res = {
            '[I,H]': H(Iop(f)) - Iop(H(f)),
            "[I',H]": H(Ipop(f)) - Ipop(H(f)),
            '[K,K+]-K+': K(Kp(f)) - Kp(K(f)) - Kp(f),
            '[K,K-]+K-': K(Km(f)) - Km(K(f)) + Km(f),
            'K+K- - Phi(H,K)': Kp(Km(f)) - Phi(0)(f),
            'K-K+ - Phi(H,K+1)': Km(Kp(f)) - Phi(1)(f),
        }
        bad = [k for k, v in res.items() if sp.expand(v) != 0]
        if bad:
            ok = False
            print(f"N={N} f={f}: FAIL {bad}")
    # form (2.6): normal-ordered coefficients of I - p2^2 and I' - p1^2 (p = -i d)
    p2sq = lambda f: -sp.diff(f, q2, 2)
    p1sq = lambda f: -sp.diff(f, q1, 2)
    form_ok = True
    for name, L in (("I - p2^2", lambda f: Iop(f) - p2sq(f)), ("I' - p1^2", lambda f: Ipop(f) - p1sq(f))):
        c = coefficients(L, 2 * N + 2)
        for (a, b), cv in c.items():
            if cv == 0:
                continue
            poly = sp.Poly(cv, q1, q2)
            degs = {sum(m) for m in poly.monoms()}
            # grading zero: coefficient degree equals derivative order; degree 1..2N: order 1..N
            if degs != {a + b} or a + b == 0 or a + b > N:
                form_ok = False
                print(f"N={N} {name}: coefficient of d1^{a} d2^{b} has degrees {degs}")
    print(f"N={N}: identities on {len(tests)} test monomials {'OK' if not bad else 'FAIL'}; form (2.6) {'OK' if form_ok else 'FAIL'}")
    ok &= form_ok
# Source reproduction at N = 2: the source's (3.1) symmetries and (3.4) ladders satisfy (3.5) = (5.1).
g1, g2 = sp.symbols('g1 g2')
H2 = add(comp(P1, P1), comp(P2, P2), poly_op([0, g1, g2], D))
r2 = lambda f: Q1(Q1(f)) + Q2(Q2(f))
I2 = lambda f: P2(P2(f)) + g1 * Q2(P2(f)) + g2 * (r2(P2(P2(f))) - C(C(f)))
I2p = lambda f: P1(P1(f)) + g1 * Q1(P1(f)) + g2 * r2(P1(P1(f)))
K1 = C
K2 = lambda f: (I2p(f) - I2(f)) / 2
K3 = lambda f: K1(K2(f)) - K2(K1(f))
Ks = lambda f: K1(f) / 2
Kp2 = lambda f: K2(f) + K3(f) / 2 - g2 * K1(K1(f)) / 2
Km2 = lambda f: K2(f) - K3(f) / 2 - g2 * K1(K1(f)) / 2
Phi1 = lambda f: (H2(f) - (2 * I) * g1 * Ks(f) - (2 * I) ** 2 * g2 * Ks(Ks(f))) / 4
KmId = lambda f: Ks(f) - f
Phi2 = lambda f: H2(f) - (-2 * I) * g1 * KmId(f) - (-2 * I) ** 2 * g2 * KmId(KmId(f))
src_ok = True
for f in tests:
    for name, v in (('[I2,H2]', H2(I2(f)) - I2(H2(f))), ("[I2',H2]", H2(I2p(f)) - I2p(H2(f))),
                    ('H2 - I2 - I2\'', H2(f) - I2(f) - I2p(f)),
                    ('K+K- - Phi1 Phi2', Kp2(Km2(f)) - Phi1(Phi2(f)))):
        if sp.expand(v) != 0:
            src_ok = False
            print(f"source N=2 f={f}: FAIL {name}")
print(f"source N=2: (3.1) symmetries, (3.2), and (3.4) ladders with (3.5)=(5.1) {'reproduced' if src_ok else 'FAIL'}")
ok &= src_ok
print("ALL_CHECKS", ok)
raise SystemExit(0 if ok else 1)
