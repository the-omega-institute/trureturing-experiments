# mpmath check (numerical evidence only): closed forms vs besseli*besselk, and the sign of
# P(-1/2)P(3/2) - P(1/2)^2 for u in {0.1, 0.5, 1, 2, 3, 5, 10}; midpoint log-convexity on several order pairs.
from mpmath import mp, mpf, besseli, besselk, exp, cosh, sinh, log
mp.dps = 40
P = lambda n, u: besseli(n, u) * besselk(n, u)
for u in [mpf(v) for v in ('0.1', '0.5', '1', '2', '3', '5', '10')]:
    Pm = (1 + exp(-2*u)) / (2*u); P0 = (1 - exp(-2*u)) / (2*u); P1 = exp(-u) * (cosh(u) - sinh(u)/u) * (1 + 1/u) / u
    err = abs(Pm - P(mpf(-0.5), u)) + abs(P0 - P(mpf(0.5), u)) + abs(P1 - P(mpf(1.5), u))
    print('u', float(u), 'closed-form err', mp.nstr(err, 3), 'P(-1/2)P(3/2)-P(1/2)^2 =', mp.nstr(Pm*P1 - P0**2, 8))
for u in [mpf(v) for v in ('0.001', '0.01', '0.1', '0.3', '1')]:
    row = []
    for (n1, n3) in [('-0.99', '-0.9'), ('-0.98', '-0.8'), ('-0.95', '-0.5'), ('-0.9', '0'), ('-0.5', '0.5'), ('0', '1')]:
        n1, n3 = mpf(n1), mpf(n3); n2 = (n1 + n3) / 2
        row.append(mp.nstr(log(P(n1, u)) + log(P(n3, u)) - 2*log(P(n2, u)), 4))
    print('u', float(u), 'midpoint log-convexity gaps', row)
