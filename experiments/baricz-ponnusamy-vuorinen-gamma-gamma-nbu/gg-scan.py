# mpmath check (numerical evidence only) for the gamma-gamma NBU question, a = b = 1, alpha = 1:
# total mass of the source density, the closed form S(u) = 2 sqrt(u) K_1(2 sqrt u), and S(2u) - S(u)^2 > 0 on a grid.
from mpmath import mp, mpf, besselk, sqrt, quad, inf, gamma
mp.dps = 30
a = b = mpf(1); al = mpf(1)
f = lambda u: 2 * (a*b/al)**((a+b)/2) * u**((a+b)/2 - 1) * besselk(a - b, 2*sqrt(a*b*u/al)) / (gamma(a)*gamma(b))
print('total mass', quad(f, [0, 1, inf]))
Sv = lambda u: 2*sqrt(u)*besselk(1, 2*sqrt(u))
print('S(1/2) via density', quad(f, [mpf('0.5'), mpf('1.5'), inf]), 'closed form', Sv(mpf('0.5')))
for u in [mpf(v) for v in ('0.01', '0.1', '0.5', '1', '2', '5')]:
    print('u', float(u), 'S(2u) - S(u)^2 =', mp.nstr(Sv(2*u) - Sv(u)**2, 6))
