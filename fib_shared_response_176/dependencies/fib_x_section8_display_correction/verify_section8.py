"""Exact outward certificates for FIB Self Calibrating X, section 8.

Inputs are the displayed rational data and VI (8.1), (8.2), (20.11),
(24.2)--(24.5). Python >=3.12 and the already installed SymPy are required.
All interval endpoints are integers divided by 10**80; every operation rounds
outward using integer division. No floating point arithmetic is used.
"""

import argparse
import ast
import hashlib
import json
import re
import sys
from fractions import Fraction as F
from math import comb
from pathlib import Path

import sympy as sp

assert sys.version_info >= (3, 12)
SCALE = 10**80
THETA, XI = F(3, 7), F(-35, 29)


def ceildiv(a, b):
    return -((-a) // b)


class Box:
    def __init__(self, lo, hi=None, raw=False):
        if isinstance(lo, Box):
            self.lo, self.hi = lo.lo, lo.hi
        elif raw:
            self.lo, self.hi = lo, lo if hi is None else hi
        else:
            lo, hi = F(lo), F(lo if hi is None else hi)
            self.lo = lo.numerator * SCALE // lo.denominator
            self.hi = ceildiv(hi.numerator * SCALE, hi.denominator)
        assert self.lo <= self.hi

    def __add__(self, other):
        other = Box(other)
        return Box(self.lo + other.lo, self.hi + other.hi, raw=True)

    __radd__ = __add__

    def __neg__(self):
        return Box(-self.hi, -self.lo, raw=True)

    def __sub__(self, other):
        return self + -Box(other)

    def __rsub__(self, other):
        return Box(other) + -self

    def __mul__(self, other):
        other = Box(other)
        products = [a * b for a in (self.lo, self.hi) for b in (other.lo, other.hi)]
        return Box(min(products) // SCALE, ceildiv(max(products), SCALE), raw=True)

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = Box(other)
        assert other.lo > 0 or other.hi < 0, "interval divisor contains zero"
        return self * Box(SCALE**2 // other.hi, ceildiv(SCALE**2, other.lo), raw=True)

    def __rtruediv__(self, other):
        return Box(other) / self

    def __pow__(self, power):
        assert isinstance(power, int) and power >= 0
        if power == 0:
            return Box(1)
        values = [self.lo**power, self.hi**power]
        lower = 0 if power % 2 == 0 and self.lo <= 0 <= self.hi else min(values)
        return Box(lower // SCALE**(power - 1),
                   ceildiv(max(values), SCALE**(power - 1)), raw=True)

    def absmax(self):
        return F(max(abs(self.lo), abs(self.hi)), SCALE)

    def sign(self):
        return "+" if self.lo > 0 else "-" if self.hi < 0 else "0"

    def endpoints(self):
        return [F(self.lo, SCALE), F(self.hi, SCALE)]

    def outward(self, digits=12):
        scale = 10**digits
        return [str(F(self.lo * scale // SCALE, scale)),
                str(F(ceildiv(self.hi * scale, SCALE), scale))]


class Dual:
    def __init__(self, value, gradient=None):
        self.value = value
        self.gradient = [0 * value] * 4 if gradient is None else gradient

    @staticmethod
    def lift(value):
        return value if isinstance(value, Dual) else Dual(value)

    def __add__(self, other):
        other = Dual.lift(other)
        return Dual(self.value + other.value,
                    [a + b for a, b in zip(self.gradient, other.gradient)])

    __radd__ = __add__

    def __neg__(self):
        return Dual(-self.value, [-a for a in self.gradient])

    def __sub__(self, other):
        return self + -Dual.lift(other)

    def __rsub__(self, other):
        return Dual.lift(other) + -self

    def __mul__(self, other):
        other = Dual.lift(other)
        return Dual(self.value * other.value,
                    [a * other.value + self.value * b
                     for a, b in zip(self.gradient, other.gradient)])

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = Dual.lift(other)
        return Dual(self.value / other.value,
                    [(a * other.value - self.value * b) / other.value**2
                     for a, b in zip(self.gradient, other.gradient)])

    def __rtruediv__(self, other):
        return Dual.lift(other) / self

    def __pow__(self, power):
        if power == 0:
            return Dual(1 + 0 * self.value)
        return Dual(self.value**power,
                    [power * self.value**(power - 1) * a for a in self.gradient])


def source(s, t):
    den = 1 + s*t
    u = (s + t + s*t) / den
    n = 1 + (2*s + t) / den
    a = (1 + s + 2*t + 3*s*t) / den
    h, ell = u*(u + 1), u*u + 3*u - 4
    delta = n*h + a*n*ell - 4*a*h
    phi = (3*h - 10)*a + (8*u*u + 18*u - 30)*n - 23*h
    v, x = 10*u*u + 19*u - 36, 20*u*u + 34*u - 70
    k, q = a/n, 5*delta/(n*phi)
    mu = ((n + 2*a)*v - 10*a*ell - 29*h)/phi
    eta = (4*(2*n + 3*a)*x - 105*a*ell - 381*h)/phi
    b0 = 120*h/phi
    gamma = b0*(1 + (1 + 2*n)*(v/h - 20/n))
    zeta = (b0*(1 + (n + a)*(x/h - 35/n)) - 210*q)/30
    d = XI*b0
    g7 = 210/n + d*(x/h - 35/n)
    mj, jm = d + (4*n - u + 1)*g7, d + (u + 2*a)*g7
    cubic = ((6 - d)/n + 3*d/(1 + u) - 3*d/u + d)/6
    return locals()


def equations(coords, arm):
    p, q = source(*coords[:2]), source(*coords[2:])
    first = [z['k'] + XI*z['q'] - THETA for z in (p, q)]
    if arm == 'M':
        return first + [p['mu'] - q['mu'], p['eta'] - q['eta']]
    return first + [1/p['n'] - 1/q['n'] + XI*(p['gamma'] - q['gamma'])/120,
                    p['zeta'] - q['zeta']]


def inverse(matrix):
    size = len(matrix)
    rows = [list(row) + [F(i == j) for j in range(size)]
            for i, row in enumerate(matrix)]
    for j in range(size):
        pivot = next(i for i in range(j, size) if rows[i][j])
        rows[j], rows[pivot] = rows[pivot], rows[j]
        val = rows[j][j]
        rows[j] = [a/val for a in rows[j]]
        for i in range(size):
            if i != j:
                val = rows[i][j]
                rows[i] = [a - val*b for a, b in zip(rows[i], rows[j])]
    return [row[size:] for row in rows]


def group(text, start):
    assert text[start] == '{'
    depth = 1
    for end in range(start + 1, len(text)):
        depth += (text[end] == '{') - (text[end] == '}')
        if depth == 0:
            return text[start + 1:end], end + 1
    raise ValueError("unclosed group")


def rational(text):
    while match := re.search(r'\\(?:d)?frac', text):
        numerator, end = group(text, match.end())
        denominator, end = group(text, end)
        text = (text[:match.start()] + '(' + str(rational(numerator)) + ')/(' +
                str(rational(denominator)) + ')' + text[end:])
    text = text.replace(r'\cdot', '*').replace('{', '(').replace('}', ')').replace('^', '**')

    def walk(node):
        if isinstance(node, ast.Constant) and type(node.value) is int:
            return F(node.value)
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
            return -walk(node.operand)
        assert isinstance(node, ast.BinOp), ast.dump(node)
        a, b = walk(node.left), walk(node.right)
        if isinstance(node.op, ast.Div):
            return a/b
        if isinstance(node.op, ast.Mult):
            return a*b
        if isinstance(node.op, ast.Pow):
            assert b.denominator == 1
            return a**b.numerator
        raise ValueError(ast.dump(node))

    return walk(ast.parse(text.strip(), mode='eval').body)


def display(section, tag):
    return section[:section.index(r'\tag{' + tag + '}')].rsplit('$$', 1)[1]


def intervals(text):
    text = text.replace(r'\left', '').replace(r'\right', '')
    return [(rational(a), rational(b))
            for a, b in re.findall(r'\[([^\[\],]+),([^\[\],]+)\]', text)]


def fractions(text):
    vals = []
    for match in re.finditer(r'\\(?:d)?frac', text):
        num, end = group(text, match.end())
        den, end = group(text, end)
        vals.append(rational(num) / rational(den))
    return vals


def polybox(poly, variables, ranges):
    """Nested Horner, using exact rational polynomial coefficients."""
    if not variables:
        return Box(F(poly))
    var, *rest = variables
    p = sp.Poly(poly, var)
    out = Box(0)
    for coefficient in p.all_coeffs():
        out = out * ranges[var] + polybox(coefficient, rest, ranges)
    return out


def centered_poly(poly, variables, ranges):
    """Exact translation before Horner limits dependency overestimation."""
    centers = {v: sum(ranges[v].endpoints())/2 for v in variables}
    translated = sp.Poly(poly, *variables, domain=sp.QQ).shift_list(
        [sp.Rational(centers[v]) for v in variables])
    return polybox(translated.as_expr(), variables,
                   {v: ranges[v] - centers[v] for v in variables})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('document', type=Path)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    data = args.document.read_bytes()
    section = data[68414:].decode()
    assert section.startswith('\n## 8. ')
    checks, measurements = [], {}

    def check(name, ok, **detail):
        checks.append({'name': name, 'passed': bool(ok), **detail})
        if not ok:
            print('FAIL', name, json.dumps(detail), flush=True)

    def encloses(name, shown, actual):
        low, high = actual.endpoints()
        check(name, shown[0] < low <= high < shown[1],
              displayed=list(map(str, shown)), computed=actual.outward(18))

    print('checking anchors', flush=True)
    centers = re.findall(r'\d+\.\d+', display(section, '8.6'))
    assert len(centers) == 8
    constants = fractions(display(section, '8.9'))
    assert len(constants) == 4
    phi_boxes = intervals(display(section, '8.10'))
    cubic_boxes = intervals(display(section, '8.11'))
    diffs = fractions(display(section, '8.12'))
    b_ends = intervals(display(section, '8.13'))[0]
    check('8.20 repeats the same selected b interval',
          tuple(fractions(display(section, '8.20'))) == b_ends)
    tight_mu = fractions(section[section.index('证明。定理8.3的 MM 压缩盒中'):])[:4]
    check('8.5 parent clock values',
          24*THETA == F(72, 7) and 60*(1 + THETA + XI) == F(2700, 203))
    for arm_index, arm in enumerate(('M', 'J')):
        c = list(map(F, centers[4*arm_index:4*arm_index + 4]))
        radius = F(1, 10**11)
        exact = equations([Dual(x, [F(i == j) for j in range(4)])
                           for i, x in enumerate(c)], arm)
        jac = [x.gradient for x in exact]
        preconditioner = inverse(jac)
        check(f'8.9 {arm} exact inverse', all(
            sum(preconditioner[i][k]*jac[k][j] for k in range(4)) == F(i == j)
            for i in range(4) for j in range(4)))
        beta = max(abs(sum(preconditioner[i][k]*exact[k].value for k in range(4)))
                   for i in range(4))
        check(f'8.9 {arm} beta', beta < constants[2*arm_index])
        coordinate_boxes = [Box(x - radius, x + radius) for x in c]
        jets = equations([Dual(x, [Box(int(i == j)) for j in range(4)])
                          for i, x in enumerate(coordinate_boxes)], arm)
        alpha = max(sum((Box(int(i == j)) - sum(
            (preconditioner[i][k]*jets[k].gradient[j] for k in range(4)), Box(0))).absmax()
                        for j in range(4)) for i in range(4))
        check(f'8.9 {arm} alpha', alpha < constants[2*arm_index + 1])
        check(f'8.9 {arm} self-map and contraction', beta + alpha*radius < radius and alpha < 1)
        check(f'8.9 {arm} root within 1e-20', beta/(1 - alpha) < F(1, 10**20))
        check(f'8.6 {arm} positive distinct sources',
              min(c) > radius and c[0] - radius > c[2] + radius)
        zs = [source(*coordinate_boxes[:2]), source(*coordinate_boxes[2:])]
        for i, z in enumerate(zs):
            encloses(f'8.10 {arm}{i + 1} phi', phi_boxes[2*arm_index + i], z['phi'])
            encloses(f'8.11 {arm}{i + 1} cubic', cubic_boxes[2*arm_index + i], z['cubic'])
            check(f'{arm}{i + 1} denominator legality', all(z[k].sign() == '+'
                  for k in ('den', 'u', 'n', 'a', 'h')) and z['phi'].sign() == '-')
        terminal = 'mj' if arm == 'M' else 'jm'
        encloses(f'8.12 {arm} alternative difference',
                 (-diffs[2*arm_index], -diffs[2*arm_index + 1]), zs[0][terminal] - zs[1][terminal])
        if arm == 'M':
            for i in range(2):
                tiny = [Box(x - F(1, 10**20), x + F(1, 10**20)) for x in c[2*i:2*i + 2]]
                mu = source(*tiny)['mu']
                encloses(f'8.20 tight mu{i + 1}', tuple(tight_mu[2*i:2*i + 2]), mu)
                encloses(f'8.13 selected b from mu{i + 1}', b_ends, mu)
        measurements[arm] = {'beta': Box(beta).outward(30), 'alpha': Box(alpha).outward(15),
                             'root_distance_bound': Box(beta/(1 - alpha)).outward(30)}

    print('building elimination polynomials from VI matrices', flush=True)
    u, b, n, a = sp.symbols('u b n a')
    h, ell = u*(u + 1), u*u + 3*u - 4
    al, be = 3*h - 10, 8*u*u + 18*u - 30
    phi = al*a + be*n - 23*h
    delta = n*h + a*n*ell - 4*a*h
    A = (10 - 3*b)*u*u + (8 - 3*b)*u + 10*b - 32
    B = (10 - 8*b)*u*u + (19 - 18*b)*u + 30*b - 36
    C = (23*b - 29)*h
    Q = sp.Matrix([[1, -sp.Rational(1, 2), -u - 1],
                   [-sp.Rational(1, 2), -1, 2*u + 1], [-u - 1, 2*u + 1, -u*u + u + 1]])
    T = sp.Matrix([[-al, (THETA*al - be - 5*XI*ell)/2, (23 + 20*XI)*h/2],
                   [(THETA*al - be - 5*XI*ell)/2, THETA*be, -(23*THETA + 5*XI)*h/2],
                   [(23 + 20*XI)*h/2, -(23*THETA + 5*XI)*h/2, 0]])
    vec = sp.Matrix([A, B, C])
    rho = sp.expand((vec.T*Q.adjugate()*vec)[0])
    sigma = sp.expand((vec.T*T.adjugate()*vec)[0])
    tau = sp.expand((vec.T*((Q + T).adjugate() - Q.adjugate() - T.adjugate())*vec)[0])
    R = sp.Poly(tau*tau - 4*rho*sigma, u, b)
    p = sp.Matrix([[n + a - 3*u - 1, n + 2*a - 4*u - 3],
                   [1 + 2*u - a, 2 + 3*u - n - a]])
    x = sp.Matrix([a, n, 1])
    check('VI 24.3 source determinant', sp.expand((x.T*Q*x)[0] - p.det()) == 0)
    check('VI 24.4 parent polynomial', sp.expand((x.T*T*x)[0] - ((THETA*n - a)*phi - 5*XI*delta)) == 0)
    psi6 = (n + 2*a)*(10*u*u + 19*u - 36) - 10*a*ell - 29*h
    check('VI 24.2 actual M response', sp.expand(A*a + B*n + C - (psi6 - b*phi)) == 0)
    # Independent literal action matrices and four-node interpolation establish
    # the formulas for the same parent history and the actual sixth action.
    M, J = sp.Matrix([[0, 1], [1, 1]]), sp.Matrix([[0, 1], [1, 0]])
    acts = [sp.eye(2)]
    for action in (M, M, J, M, M):
        acts.append(action*acts[-1])
    trace = lambda action: sp.expand(sp.trace(action*p))
    check('8.1 source traces on literal MMJMM parent',
          [trace(action) for action in acts] == [1, u, u + 1, n, a, n + a])
    check('8.3/8.4 actual sixth and seventh trace factors',
          [trace(M*acts[-1]), trace(M*M*acts[-1]), trace(J*M*acts[-1]),
           trace(J*acts[-1]), trace(J*J*acts[-1]), trace(M*J*acts[-1])] ==
          [n + 2*a, 2*n + 3*a, 4*n - u + 1, 1 + 2*n, n + a, u + 2*a])
    drift = sp.Symbol('drift')
    nodes = [-drift, -drift/u, -drift/(u + 1), (6 - drift)/n]

    def interpolate(clock):
        return sum(value*sp.prod(sp.Rational(clock - j, i - j)
                                for j in range(4) if j != i)
                   for i, value in enumerate(nodes))

    V6, X7 = 10*u*u + 19*u - 36, 20*u*u + 34*u - 70
    g6, g7 = 120/n + drift*(V6/h - 20/n), 210/n + drift*(X7/h - 35/n)
    check('8.2 gain from whole cubic interpolation',
          sp.cancel(interpolate(6) - g6) == 0 and sp.cancel(interpolate(7) - g7) == 0)
    theta_local, xi_local = a/n + drift*delta/(24*h*n), drift*phi/(120*h)
    q_src = 5*delta/(n*phi)
    psi7 = 4*(2*n + 3*a)*X7 - 105*a*ell - 381*h
    gamma = 120*h/phi*(1 + (1 + 2*n)*(V6/h - 20/n))
    zeta = (120*h/phi*(1 + (n + a)*(X7/h - 35/n)) - 210*q_src)/30
    identities = [
        (n + 2*a)*g6 + drift - (120*(1 + 2*theta_local) + 120*xi_local*psi6/phi),
        (2*n + 3*a)*g7 + drift - (210*(2 + 3*theta_local) + 30*xi_local*psi7/phi),
        (1 + 2*n)*g6 + drift - (240 + 120/n + xi_local*gamma),
        (n + a)*g7 + drift - (210*(1 + theta_local) + 30*xi_local*zeta),
        a*interpolate(4) + drift - 24*theta_local,
        (n + a)*interpolate(5) + drift - 60*(1 + theta_local + xi_local),
    ]
    check('8.1--8.7 response and parent identities on full source coordinates',
          all(sp.cancel(identity) == 0 for identity in identities))
    xa = sp.Matrix([-B*n - C, A*n, A])
    qa = sp.Poly((xa.T*Q*xa)[0], n)
    ta = sp.Poly((xa.T*T*xa)[0], n)
    q2, q1, q0 = qa.all_coeffs()
    t2, t1, t0 = ta.all_coeffs()
    N, D = sp.expand(t2*q0 - q2*t0), sp.expand(q2*t1 - t2*q1)
    nhat, dhat = [sp.cancel(-203*z/A) for z in (N, D)]
    check('8.17 polynomial factors and degrees',
          sp.denom(nhat) == 1 and sp.denom(dhat) == 1 and
          (sp.degree(nhat, u), sp.degree(dhat, u), sp.degree(nhat, b), sp.degree(dhat, b)) == (10, 9, 3, 3))
    chart_resultant = sp.expand((q2*t0 - q0*t2)**2 - (q2*t1 - q1*t2)*(q1*t0 - q0*t1))
    check('8.17 chart resultant equals A^4 R', sp.Poly(chart_resultant - A**4*R.as_expr(), u, b).is_zero)
    check('8.17 reconstructed n solves both quadratics at every R root', all(
        sp.div(sp.Poly(lead*nhat**2 + middle*nhat*dhat + const*dhat**2, u, b), R)[1].is_zero
        for lead, middle, const in ((q2, q1, q0), (t2, t1, t0))))
    check('8.16 degree 16', R.degree(u) == 16)
    b_box = Box(*b_ends)

    print('checking all 13 Descartes rows, including the infinite tail', flush=True)
    boundaries = [F(v) for v in ['0', '1821/5000', '3643/10000', '5607/5000',
                  '2243/2000', '11603/10000', '2901/2500', '1177/1000', '11771/10000',
                  '3199/2000', '15997/10000', '41691/5000', '16677/2000']]
    expected = ['+'*17, '+'*3 + '-'*14, '-'*17, '-'*3 + '+'*14, '+'*17,
                '+'*10 + '-'*7, '-'*17, '-'*7 + '+'*10, '+'*17,
                '+'*7 + '-'*10, '-'*17, '-'*7 + '+'*10, '+'*17]
    table16 = display(section, '8.16')
    shown_signs = re.findall(r'[\])]&([^&]+)&', table16)
    expanded = [''.join(s*int(k or 1) for s, k in re.findall(r'([+-])(?:\^\{?(\d+)\}?)?', word))
                for word in shown_signs]
    check('8.16 displayed sign strings', expanded == expected)
    displayed_bounds = re.findall(r'\[([^,\]]+),([^\])]+)[\])]', table16)
    check('8.16 displayed interval partition', len(displayed_bounds) == 13 and all(
        rational(left) == boundaries[i] and
        (right == r'\infty' if i == 12 else rational(right) == boundaries[i + 1])
        for i, (left, right) in enumerate(displayed_bounds)))
    coeff = [[F(R.coeff_monomial(u**k*b**j)) for k in range(17)] for j in range(R.degree(b) + 1)]
    for row, lo in enumerate(boundaries):
        transformed = []
        for power in range(17):
            by_b = []
            for j, coefficients in enumerate(coeff):
                val = F(0)
                for k, c in enumerate(coefficients):
                    if row < 12:
                        hi = boundaries[row + 1]
                        val += c*sum(F(comb(k, e)*comb(16 - k, power - e))*lo**(k - e)*hi**e
                                     for e in range(max(0, power - 16 + k), min(k, power) + 1))
                    elif k >= 16 - power:
                        val += c*comb(k, 16 - power)*lo**(k - 16 + power)
                by_b.append(val)
            val = Box(0)
            for c in reversed(by_b):
                val = val*b_box + c
            transformed.append(val.sign())
        signs = ''.join(transformed)
        check(f'8.16 row {row + 1} all 17 coefficients', signs == expected[row], observed=signs)
        endpoint_signs = [polybox(R.as_expr().subs(u, sp.Rational(lo)), [b], {b: b_box}).sign()]
        endpoint_signs += [polybox(R.as_expr().subs(u, sp.Rational(boundaries[row + 1])), [b], {b: b_box}).sign()
                           if row < 12 else polybox(sp.Poly(R.as_expr(), u).LC(), [b], {b: b_box}).sign()]
        check(f'8.16 row {row + 1} endpoint signs', endpoint_signs == [expected[row][0], expected[row][-1]])

    print('checking root boxes, all reconstructed sources and terminal responses', flush=True)
    cells = intervals(display(section, '8.18'))
    assert len(cells) == 30
    shown_source_signs = re.findall(r'&([+-]{4})(?:\\\\|\s*\n\\end)', display(section, '8.18'))
    assert len(shown_source_signs) == 6
    terminal_boxes = intervals(display(section, '8.15'))
    root_results = []
    for i in range(6):
        shown = cells[5*i:5*i + 5]
        ubox = Box(*shown[0])
        ranges = {u: ubox, b: b_box}
        root_end_signs = [polybox(R.as_expr().subs(u, sp.Rational(ep)), [b], {b: b_box}).sign()
                          for ep in shown[0]]
        check(f'8.18 root {i + 1} inside unique broad root box',
              boundaries[2*i + 1] < shown[0][0] < shown[0][1] < boundaries[2*i + 2] and
              root_end_signs == [expected[2*i + 1][0], expected[2*i + 1][-1]])
        av = centered_poly(A, [u, b], ranges)
        bv = centered_poly(B, [u, b], ranges)
        cv = centered_poly(C, [u, b], ranges)
        dv = centered_poly(dhat, [u, b], ranges)
        nv = centered_poly(nhat, [u, b], ranges)/dv
        a_val = -(bv*nv + cv)/av
        for label, interval, val in zip(('A', 'Dhat', 'n', 'a'), shown[1:], (av, dv, nv, a_val)):
            encloses(f'8.18 root {i + 1} {label}', interval, val)
        lead_q = centered_poly(q2, [u, b], ranges)
        lead_t = centered_poly(t2, [u, b], ranges)
        original_D = -av*dv/203
        check(f'8.18 root {i + 1} pivot, leading and linear coefficients nonzero',
              all(val.sign() != '0' for val in (av, dv, lead_q, lead_t, original_D)))
        source_values = [nv + a_val - 3*ubox - 1, nv + 2*a_val - 4*ubox - 3,
                         1 + 2*ubox - a_val, 2 + 3*ubox - nv - a_val]
        signs = ''.join(val.sign() for val in source_values)
        check(f'8.18 root {i + 1} source signs', signs == shown_source_signs[i], observed=signs)
        phi_val = (3*ubox*(ubox + 1) - 10)*a_val + (8*ubox**2 + 18*ubox - 30)*nv - 23*ubox*(ubox + 1)
        check(f'8.18 root {i + 1} drift/gain denominator legality', phi_val.sign() != '0' and nv.sign() != '0')
        result = {'root': i + 1, 'source_signs': signs,
                  'A': av.outward(), 'Dhat': dv.outward(), 'n': nv.outward(), 'a': a_val.outward()}
        if i in (1, 3):
            drift = 120*ubox*(ubox + 1)*XI/phi_val
            gain7 = 210/nv + drift*((20*ubox**2 + 34*ubox - 70)/(ubox*(ubox + 1)) - 35/nv)
            mj = drift + (4*nv - ubox + 1)*gain7
            encloses(f'8.15 root {i + 1} actual M then J', terminal_boxes[0 if i == 1 else 1], mj)
            result['actual_MJ'] = mj.outward()
        root_results.append(result)
    check('8.15 terminal intervals disjoint', terminal_boxes[1][1] < terminal_boxes[0][0])
    measurements['roots'] = root_results
    source_document = args.document.with_name('FIB_SELF_CALIBRATING_RELATION_MATHEMATICS_VI.md')
    report = {'document': str(args.document.resolve()), 'document_sha256': hashlib.sha256(data).hexdigest(),
              'source_document': str(source_document.resolve()),
              'source_document_sha256': hashlib.sha256(source_document.read_bytes()).hexdigest(),
              'verifier_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'python': sys.version.split()[0], 'sympy': sp.__version__, 'interval_grid': str(SCALE),
              'source_formulas': ['VI (8.1), (8.2), (8.5)--(8.9)', 'VI (13.1)--(13.3)',
                                  'VI (20.11)', 'VI (24.2)--(24.5)'],
              'arithmetic': 'rational symbolic algebra and outward integer interval arithmetic; no floats',
              'checks': checks, 'measurements': measurements,
              'passed': all(c['passed'] for c in checks),
              'not_established': ['Lean verification', 'general common-parent L2',
                                  'other sixth-response fibers', 'root eight/nine-read decision or nine-read upper bound']}
    if args.output:
        args.output.write_text(json.dumps(report, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps({'passed': report['passed'], 'checks': len(checks),
                      'failures': [c for c in checks if not c['passed']], 'measurements': measurements}, indent=2), flush=True)
    return 0 if report['passed'] else 1


if __name__ == '__main__':
    sys.exit(main())
