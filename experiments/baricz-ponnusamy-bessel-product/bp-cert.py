# Exact rational certificate (standard library only).
# Baricz-Ponnusamy, Proc. AMS 141 (2013), arXiv:1010.3346, Conjecture: nu -> I_nu(u) K_nu(u) is log-convex on (-1, oo).
# Log-convexity forces P(-1/2) P(3/2) >= P(1/2)^2 for P(nu) = I_nu(u) K_nu(u).  At half-integer orders
#   I_{-1/2}K_{-1/2}(u) = (1 + e^{-2u})/(2u),  I_{1/2}K_{1/2}(u) = (1 - e^{-2u})/(2u),
#   I_{3/2}K_{3/2}(u) = e^{-u}(cosh u - sinh u/u)(1 + 1/u)/u,
# so at u = 1 (cosh 1 - sinh 1 = e^{-1}): P(-1/2) = (1 + t)/2, P(1/2) = (1 - t)/2, P(3/2) = 2t with t = e^{-2}.
# The inequality to refute is (1 + t) t >= (1 - t)^2 / 4.  The difference D(t) = (1+t)t - (1-t)^2/4 is
# increasing on [0, 1] (D'(t) = 1 + 2t + (1-t)/2 > 0), so it suffices to bound t = e^{-2} from above.
from fractions import Fraction as F
e_lo = F(2718281828, 10**9)           # e > 2.718281828  (e = 2.71828182845...)
t_hi = 1 / e_lo**2                    # t = e^{-2} < t_hi
D = lambda t: (1 + t) * t - (1 - t)**2 / 4
assert D(t_hi) < 0
print("t = e^-2 <", float(t_hi), "; D(t) <= D(t_hi) =", float(D(t_hi)), "< 0  (exact)")
