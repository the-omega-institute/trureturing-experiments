#!/usr/bin/env python3
"""Check: the triple sum S of Cork and Disney-Hogg (arXiv:2505.00075v2, Remark 13, Eq. (18)) is positive for all
nonzero, not-all-equal y_i, through the identity S = A * sum_i q_i |y_i - m|^2 + |A m - B|^2 / t.
Part A (sympy, exact rationals): the identity on random rational instances in R^2, R^3 and R^4.
Part B (floating point): positivity and the identity on 20 000 random configurations, the equality case.
Exit 0 iff every check passes."""
import random, sys
import sympy as sp
fails = []
def check(name, ok):
    print(('PASS ' if ok else 'FAIL ') + name)
    if not ok: fails.append(name)
def dot(u, v): return sum(a * b for a, b in zip(u, v))
def sub(u, v): return [a - b for a, b in zip(u, v)]
def S_triple(lam, mu, ys):
    n = len(ys); w = [l * l for l in lam]; r = [dot(y, y) for y in ys]; t = mu * mu
    rho0 = sum(w[i] / r[i] for i in range(n))
    tot = 0
    for i in range(n):
        for j in range(n):
            for k in range(n):
                tot += w[i] * w[j] * w[k] * dot(sub(ys[j], ys[i]), sub(ys[j], ys[k])) / ((r[i] + t) * (r[j] + t) * r[j] * r[k])
    return tot / rho0
def S_sos(lam, mu, ys):
    n = len(ys); dim = len(ys[0]); w = [l * l for l in lam]; r = [dot(y, y) for y in ys]; t = mu * mu
    A = sum(w[i] / (r[i] + t) for i in range(n)); C = sum(w[i] / r[i] for i in range(n))
    B = [sum(w[i] * ys[i][c] / (r[i] + t) for i in range(n)) for c in range(dim)]
    Dv = [sum(w[i] * ys[i][c] / r[i] for i in range(n)) for c in range(dim)]
    m = [x / C for x in Dv]; q = [w[i] / (r[i] * (r[i] + t)) for i in range(n)]
    closed = A * A + (C * dot(B, B) - A * dot(Dv, Dv)) / (C * t)
    sos = A * sum(q[i] * dot(sub(ys[i], m), sub(ys[i], m)) for i in range(n)) + dot(sub([A * x for x in m], B), sub([A * x for x in m], B)) / t
    return closed, sos
rnd = random.Random(20261019)
def rq(): return sp.Rational(rnd.randint(-9, 9), rnd.randint(1, 5))
ok_exact = True
for trial in range(30):
    dim = rnd.choice([2, 3, 4]); n = rnd.randint(2, 4)
    ys = []
    while len(ys) < n:
        y = [rq() for _ in range(dim)]
        if any(y) and y not in ys: ys.append(y)
    lam = [sp.Rational(rnd.randint(1, 6), rnd.randint(1, 4)) for _ in range(n)]; mu = sp.Rational(rnd.randint(1, 6), rnd.randint(1, 4))
    St = sp.nsimplify(S_triple(lam, mu, ys)); cl, so = S_sos(lam, mu, ys)
    if sp.simplify(St - cl) != 0 or sp.simplify(St - so) != 0 or not (St > 0): ok_exact = False; print('INFO bad', trial, St, cl, so)
check('exact: S (Eq. (18)) = closed form = sum-of-squares form, and S > 0, on 30 rational instances (dim 2-4, N+1 = 2-4)', ok_exact)
worst = 1e9; ident = 0.0
for trial in range(20000):
    dim = 3; n = rnd.randint(2, 7)
    ys = [[rnd.gauss(0, rnd.choice([0.1, 1, 10])) for _ in range(dim)] for _ in range(n)]
    lam = [rnd.uniform(0.05, 5) for _ in range(n)]; mu = rnd.uniform(0.05, 5)
    St = S_triple(lam, mu, ys); cl, so = S_sos(lam, mu, ys)
    worst = min(worst, St / max(1e-300, abs(cl) + 1e-300)); ident = max(ident, abs(St - so) / max(1.0, abs(St)))
print(f'INFO float: min S/|closed| = {worst:.6f}, max relative |S - sos| = {ident:.2e}')
check('float: S > 0 on 20 000 random configurations in R^3', worst > 0)
check('float: identity holds to 1e-9 relative', ident < 1e-9)
ys = [[1, 2, 3]] * 3; lam = [1, 2, 3]
check('equality case: all y_i equal gives S = 0', abs(S_triple(lam, 1.0, ys)) < 1e-12)
print('RESULT', 'FAIL' if fails else 'PASS', len(fails))
sys.exit(1 if fails else 0)
