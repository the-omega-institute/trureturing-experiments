#!/usr/bin/env python3
"""Exploratory BA-type alternating maximization; floating point, not a certificate.

Each side is optimized with the other fixed, under its information budget.
SLSQP solves the constrained inner problems. This is not a claim to have found
the global optimum; check.py supplies separate exact finite-grid evidence.
"""
import json
import numpy as np
import scipy
from scipy.optimize import minimize

LN2 = np.log(2.0)
SEED = 20221103
rng = np.random.default_rng(SEED)


def hb(a):
    return -(a*np.log(a)+(1-a)*np.log1p(-a))/LN2


def rate(k):
    pu = k.mean(axis=0)
    return np.sum(k*np.log(k/pu))/2/LN2


def value(k, l, p):
    source = np.array([[1-p, p], [p, 1-p]])/2
    j = k.T @ source @ l
    return np.sum(j*np.log(j/(j.sum(axis=1)[:, None]*j.sum(axis=0)[None, :])))/LN2


def embed(a, n):
    left = n//2
    return np.array([[*( [(1-a)/left]*left), *([a/(n-left)]*(n-left))],
                     [*([a/left]*left), *([(1-a)/(n-left)]*(n-left))]])


def random_feasible(n, budget):
    k = rng.dirichlet(np.ones(n)*2, 2)
    center = k.mean(axis=0)
    while rate(k) > budget:
        k = (k+center)/2
    return k


def side(k, other, p, budget):
    n = k.shape[1]
    transition = np.array([[1-p, p], [p, 1-p]]) @ other
    scale = (1-2*p)**2

    def fun(z):
        return -value(z.reshape(2, n), other, p)/scale

    def jac(z):
        kk = z.reshape(2, n)
        j = kk.T @ transition / 2
        ratio = j/(j.sum(axis=1)[:, None]*j.sum(axis=0)[None, :])
        return (-((np.log(ratio) @ transition.T).T)/2/LN2/scale).ravel()

    def eq(z):
        return z.reshape(2, n).sum(axis=1)-1

    def eqjac(z):
        out = np.zeros((2, 2*n))
        out[0, :n] = out[1, n:] = 1
        return out

    def cap(z):
        return budget-rate(z.reshape(2, n))

    def capjac(z):
        kk = z.reshape(2, n)
        return (-np.log(kk/kk.mean(axis=0))/2/LN2).ravel()

    res = minimize(fun, k.ravel(), jac=jac, method='SLSQP',
                   bounds=[(1e-10, 1-1e-10)]*(2*n),
                   constraints=[{'type': 'eq', 'fun': eq, 'jac': eqjac},
                                {'type': 'ineq', 'fun': cap, 'jac': capjac}],
                   options={'ftol': 2e-13, 'maxiter': 300})
    kk = res.x.reshape(2, n)
    feasible = abs(eq(res.x)).max() <= 1e-8 and cap(res.x) >= -1e-9
    if feasible and value(kk, other, p) >= value(k, other, p)-1e-12:
        return kk, res.success
    return k, False


def main():
    print(json.dumps({'method': 'BA-type alternating constrained maximization; SLSQP inner solver',
                      'seed': SEED, 'numpy': np.__version__, 'scipy': scipy.__version__,
                      'scope': 'exploratory floating point; no global-optimum certificate'}))
    for a, b in [(1/4, 1/4), (1/8, 1/3)]:
        cu, cv = 1-hb(a), 1-hb(b)
        for p in [0.49, 0.495]:
            target = 1-hb(a*(1-p)+(1-a)*p+(1-2*(a*(1-p)+(1-a)*p))*b)
            for n in [2, 3, 4]:
                best = -1.0
                worst_rate_excess = 0.0
                inner_failures = 0
                for restart in range(6):
                    k = embed(a, n) if restart == 0 else random_feasible(n, cu)
                    l = embed(b, n) if restart == 0 else random_feasible(n, cv)
                    old = value(k, l, p)
                    for iteration in range(30):
                        k, ok1 = side(k, l, p, cu)
                        l, ok2 = side(l, k, p, cv)
                        inner_failures += int(not ok1)+int(not ok2)
                        new = value(k, l, p)
                        if abs(new-old) < 2e-14:
                            break
                        old = new
                    best = max(best, value(k, l, p))
                    worst_rate_excess = max(worst_rate_excess, rate(k)-cu, rate(l)-cv)
                print(json.dumps({'alpha': a, 'beta': b, 'Cu': cu, 'Cv': cv, 'p': p,
                                  'alphabet_sizes': [n, n], 'restarts': 6,
                                  'BSC_bits': target, 'best_bits': best,
                                  'best_minus_BSC_bits': best-target,
                                  'max_rate_excess': worst_rate_excess,
                                  'inner_unsuccessful_statuses': inner_failures}, sort_keys=True))
                assert best-target < 1e-8
                assert worst_rate_excess < 1e-8
    print('NUMERICAL_CHECKS_OK')


if __name__ == '__main__':
    main()
