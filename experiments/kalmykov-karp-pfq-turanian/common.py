"""Shared definitions: c(t) = prod Gamma(a_i+t) / prod Gamma(b_j+t) and the
coefficients delta_m of the shifted Turanian of f(mu;x) = sum_n c(mu+n) x^n/n!."""
import itertools
import mpmath as mp


def c(a, b, t):
    return mp.e ** (sum(mp.loggamma(x + t) for x in a) - sum(mp.loggamma(x + t) for x in b))


def delta(a, b, mu, al, be, m):
    return mp.fsum(
        (c(a, b, mu + al + k) * c(a, b, mu + be + m - k) - c(a, b, mu + k) * c(a, b, mu + al + be + m - k))
        / (mp.factorial(k) * mp.factorial(m - k))
        for k in range(m + 1))


def esym(v, k):
    if k < 0 or k > len(v):
        return mp.mpf(0)
    if k == 0:
        return mp.mpf(1)
    return mp.fsum(mp.fprod(s) for s in itertools.combinations(v, k))


def cond13(a, b):
    p, q = len(a), len(b)
    r = [esym(b, q - j) / esym(a, p - j) for j in range(p)] + [esym(b, q - p)]
    return all(r[i] <= r[i + 1] * (1 + mp.mpf("1e-30")) for i in range(len(r) - 1))


def f(a, b, mu, x, terms):
    return mp.fsum(c(a, b, mu + n) * x ** n / mp.factorial(n) for n in range(terms))


def turanian(a, b, mu, al, be, x, terms):
    return f(a, b, mu + al, x, terms) * f(a, b, mu + be, x, terms) - f(a, b, mu, x, terms) * f(a, b, mu + al + be, x, terms)
