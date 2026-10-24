"""Finite multiprecision grid for the function clause; Python 3 + mpmath."""
import mpmath as mp
mp.mp.dps = 80

def series(t, q, x):
    term = mp.gamma(t) ** (-q)
    value = term
    for n in range(399):
        term *= x / ((n+1) * (t+n)**q)
        value += term
    return value

for q in (9, 10, 11, 12, 13, 14, 16, 18):
    best = None
    for i in range(121):
        x = mp.power(10, mp.mpf(i)/20)
        a = series(mp.mpf(2), q, x)
        r = series(mp.mpf('1.5'), q, x)
        s = series(mp.mpf('2.5'), q, x)
        ratio = 1-r*s/a**2
        if best is None or ratio < best[0]:
            best = (ratio, x)
    print(q, mp.nstr(best[0], 12), mp.nstr(best[1], 12))
