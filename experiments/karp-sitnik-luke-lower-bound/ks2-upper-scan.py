# mpmath scan (numerical evidence only) for Karp-Sitnik Conjecture 2 with sigma in (0,1):
# a = (1, 1), b = (1/10, 2), sigma in {0.01, 0.05, 0.09}, x in {0.01, 0.1, 0.5, 2, 10, 1000}.
# Prints (q+1)Fq(sigma, a; b; -x) and the conjectured upper bound (1 + x K^(1/sigma))^(-sigma),
# K = prod Gamma(a_i)Gamma(b_i - sigma)/(Gamma(b_i)Gamma(a_i - sigma)), and whether F < bound.
from mpmath import mp, mpf, hyper, gamma
mp.dps = 30
a = [mpf(1), mpf(1)]; b = [mpf('0.1'), mpf(2)]
for s in [mpf('0.01'), mpf('0.05'), mpf('0.09')]:
    K = mpf(1)
    for ai, bi in zip(a, b):
        K *= gamma(ai) * gamma(bi - s) / (gamma(bi) * gamma(ai - s))
    for x in [mpf('0.01'), mpf('0.1'), mpf('0.5'), mpf(2), mpf(10), mpf(1000)]:
        F = hyper([s] + a, b, -x); B = 1 / (1 + x * K ** (1 / s)) ** s
        print('sigma', float(s), 'x', float(x), 'F', mp.nstr(F, 10), 'bound', mp.nstr(B, 10), 'F<bound', F < B)
