"""Exact checks of four scalar witness families and universal identities.

Symbols remain arbitrary real numbers. Raw h and g are always rebuilt from
the same D_alpha; full four-dimensional products define every budget.
The final identities support the universal inequalities stated in the data.
"""
import json
import sys

assert sys.version_info >= (3, 12)
from pathlib import Path
import sympy as s

I = s.I
alpha = s.symbols('alpha', real=True)
eps = s.symbols('epsilon', real=True)
tau = s.symbols('tau', positive=True)
M = s.symbols('M', real=True)
C = s.diag(1, 2, 3, 0)


def e(i, j):
    a = s.zeros(4)
    a[i - 1, j - 1] = 1
    return a


def F(i, j, z):
    return z * e(i, j) + s.conjugate(z) * e(j, i)


def circle(a, b):
    return (a * b + b * a) / 2


def simp(a):
    return a.applyfunc(s.simplify) if isinstance(a, s.MatrixBase) else s.simplify(a)


def eq(a, b):
    residual = simp(a - b)
    assert residual == (s.zeros(*a.shape) if isinstance(a, s.MatrixBase) else 0), residual


def K(a, b):
    return I * (a * b - b * a) / 2


def raw(x, t, v):
    h, g = alpha * K(C, x), alpha * K(C, t)
    r, k = circle(C, x), circle(C, t)
    L = -2 * circle(x, h - k) - 2 * circle(t, r + g)
    readout = simp((v.H * L * v)[0])
    rows = []
    for d in (1, -1):
        R, S = x + d * (h - k), t + d * (r + g)
        A = simp((v.H * R * R * v)[0])
        B = simp((v.H * S * S * v)[0])
        T = simp((v.H * circle(R, S) * v)[0])
        budget = s.factor(A * B - T * T)
        rows.append((d, simp(R * v), simp(S * v), budget))
    return readout, rows


def main():
    result = {}
    b = s.sqrt(8)
    v = s.Matrix([0, 0, 1, 1])
    L, rows = raw(F(1, 3, 1), F(1, 4, b) + F(2, 4, eps), v)
    eq(L, -3 * b)
    for d, Rv, Sv, budget in rows:
        A = 1 - d * (b / 2 + I * alpha)
        B = b + d * (2 + I * alpha * b / 2)
        eq(Rv, s.Matrix([A, -d * eps, 0, 0]))
        eq(Sv, s.Matrix([B, eps * (1 + d * I * alpha), 0, 0]))
        eq(s.im(s.conjugate(A) * B), d * 3 * alpha * b / 2)
        eq(budget.subs(eps, 0), 18 * alpha**2)
    result['alternative_abs_alpha_lt_one'] = {
        'exact_raw_readouts_checked': True, 'phi_L': str(L),
        'budget_at_epsilon_zero': '18*alpha**2 for each sign',
        'strict_bound': 'For |alpha|<1, epsilon=(1-|alpha|)/100. Real wedge perturbation <=17*epsilon per sign: ||dR wedge S0||<7epsilon, ||R0 wedge dS||<8epsilon, ||dR wedge dS||<2epsilon. RHS<=3sqrt(8)|alpha|+34epsilon<3sqrt(8)=LHS.'
    }
    w = 3 * (1 + I) / 2
    x = (2 - 3 * I) * e(3, 1) + (2 + 3 * I) * e(1, 3)
    x += (3 * I - 1) * e(3, 2) + (-3 * I - 1) * e(2, 3)
    v = s.Matrix([1, 1, 0, 0])
    L, rows = raw(x, F(1, 4, 2 * eps) + F(2, 4, -eps), v)
    eq(L, -3 * alpha)
    for d, Rv, Sv, budget in rows:
        eq(Rv, s.Matrix([0, 0, 1 + d * alpha * w, 0]))
        eq(Sv, s.Matrix([0, 0, d * w, eps]))
        eq(budget, s.Rational(9, 4) + eps**2 * s.Abs(1 + d * alpha * w)**2)
    result['alternative_abs_alpha_gt_one'] = {
        'exact_raw_readouts_checked': True, 'phi_L': str(L),
        'budget': '9/4+epsilon**2*|1+delta*alpha*w|**2',
        'strict_bound': 'For s=|alpha|>1, epsilon=(s-1)/(2(1+s|w|)). RHS<=3+2epsilon(1+s|w|)=s+2<3s. Both budgets >=9/4.'
    }
    v = s.Matrix([2 * I * M - 2, 2 - I * M, 0, 0])
    L, rows = raw(s.zeros(4), F(1, 4, 1) + F(2, 4, 1), v)
    eq(L, 2 * alpha * M)
    for d, Rv, Sv, budget in rows:
        eq(Rv, s.Matrix([0, 0, 0, -d]))
        eq(Sv, s.Matrix([0, 0, 0, I * (M - d * alpha)]))
        eq(budget, (M - d * alpha)**2)
    result['finite_abs_alpha_gt_one'] = {
        'exact_raw_readouts_checked': True, 'phi_L': str(L),
        'budget': '(M-delta*alpha)**2',
        'strict_bound': 'M=|alpha|+1>|alpha|. RHS=2M<2|alpha|M=LHS for |alpha|>1.'
    }
    v = s.Matrix([tau, 0, 0, 1])
    L, rows = raw(e(1, 1) / tau**2, F(1, 4, 1) + F(2, 4, tau), v)
    eq(L, -1 / tau)
    budgets = []
    for d, Rv, Sv, budget in rows:
        eq(Rv, s.Matrix([1 / tau - d / 2, -d * tau, 0, -d * tau / 2]))
        eq(Sv, s.Matrix([1 + d / tau + d * I * alpha / 2,
                            tau + d * I * alpha * tau, 0,
                            tau - d * I * alpha * tau / 2]))
        # A real-coordinate minor alone contributes 9/4 to the full real Gram.
        eq(s.re(Rv[0]) * s.re(Sv[3]) - s.re(Rv[3]) * s.re(Sv[0]), s.Rational(3, 2))
        budgets.append(budget)
    poly = s.factor(2 * sum(budgets) * tau**2)
    target = (alpha**2 + (25 + s.Rational(21, 4) * alpha**2) * tau**2
              + (1 + s.Rational(5, 2) * alpha**2) * tau**4
              + (1 + s.Rational(25, 4) * alpha**2) * tau**6)
    eq(poly, target)
    result['finite_abs_alpha_lt_one'] = {
        'exact_raw_readouts_checked': True, 'phi_L': str(L),
        'two_sum_budgets_times_tau_squared': str(s.expand(poly)),
        'strict_bound': 'For |alpha|<1 set tau^2=(1-alpha^2)/82. The polynomial C<=alpha^2+41tau^2=(1+alpha^2)/2<1. (sqrt(B+)+sqrt(B-))^2<=2(B++B-)=C/tau^2<1/tau^2=LHS^2. Both B>=9/4 by the real minor.'
    }
    ca = s.diag(1, 2, 3)
    cauchy = s.Matrix(3, 3, lambda i, j: s.Rational(2, i + j + 2))
    ones = s.ones(3, 1)
    eq(cauchy.inv() * ones, s.Matrix([6, -30, 30]))
    eq((ca * cauchy + cauchy * ca) / 2, s.ones(3))
    y, z = s.Matrix([1, -1, 0]), s.Matrix([1, 0, -1])
    d1, d2, d3 = s.symbols('d1 d2 d3', real=True)
    d = s.diag(d1, d2, d3)
    compression = s.factor((y.T * (d * cauchy**2 - cauchy**2 * d) * z)[0])
    eq(compression, s.Rational(4, 15) * (d1 - 2 * d2 + d3))
    result['cauchy_and_H_reduction'] = {
        'inverse_times_ones': [6, -30, 30],
        'determinant': str(cauchy.det()),
        'corrected_compressed_commutator': str(compression),
        'correct_test': '(e1-e2)^T[d,C^2](e1-e3)=0; gives h1-4h2+3h3=0 for d=(h1,2h2,3h3)'
    }
    # Free noncommutative-polynomial identities hold in every matrix size.
    nc_c, nc_x, nc_t = s.symbols('c x t', commutative=False)
    records = []
    for sigma in (1, -1):
        h = sigma * K(nc_c, nc_x)
        g = sigma * K(nc_c, nc_t)
        r, k = circle(nc_c, nc_x), circle(nc_c, nc_t)
        zs, products, sums, brackets = [], [], [], []
        for delta in (1, -1):
            R = nc_x + delta * (h - k)
            S = nc_t + delta * (r + g)
            Z = R + I * S
            target = ((1 + delta * I * nc_c) * (nc_x + I * nc_t)
                      if sigma == 1 else
                      (nc_x + I * nc_t) * (1 + delta * I * nc_c))
            assert s.expand(Z - target) == 0
            product = (R - I * S) * (R + I * S) if sigma == 1 else (R + I * S) * (R - I * S)
            products.append(product)
            sums.append(R * R + S * S)
            brackets.append(sigma * K(R, S))
        nc_L = -2 * circle(nc_x, h - k) - 2 * circle(nc_t, r + g)
        assert s.expand(products[0] - products[1]) == 0
        assert s.expand(sums[0] - sums[1] + 2 * nc_L) == 0
        assert s.expand(brackets[0] - brackets[1] - nc_L) == 0
        records.append(sigma)
    result['universal_sufficiency'] = {
        'checked_orientations': records,
        'method': 'Exact free noncommutative polynomial expansion; no commutativity or finite sampling assumed.',
        'correct_products': 'sigma=+1: Z_delta^* Z_delta; sigma=-1: Z_delta Z_delta^*. Each is delta independent.',
        'result': 'K_sigma(R+,S+)-K_sigma(R-,S-)=L. Positive-functional 2x2 complex Gram positivity gives |phi(K_sigma(R,S))|^2<=B_phi(R,S); triangle inequality gives raw PB for every positive complex-linear phi, including zero mass and zero budgets.'
    }
    result['status'] = 'All symbolic assertions passed; no numerical sampling or Lean verification.'
    Path(__file__).with_name('exact_matrix_checks.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
