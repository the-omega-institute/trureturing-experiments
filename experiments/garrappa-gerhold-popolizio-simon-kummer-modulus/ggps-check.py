# Witness check for Garrappa-Gerhold-Popolizio-Simon Remark 6 inequality (2.4):
# |M(a, a-b, z)| >= M(a, a-b, Re z).  M = Kummer 1F1.
import mpmath as mp
mp.mp.dps = 50
M = lambda a, b, z: mp.hyp1f1(a, b, z)
def report(tag, a, beta, z):
    b = a - beta
    lhs = abs(M(a, b, z)); rhs = M(a, b, mp.re(z))
    print(tag, "a=%s beta=%s b=%s z=%s" % (a, beta, b, mp.nstr(z, 15)),
          "|M(z)|=%s M(Re z)=%s fails=%s" % (mp.nstr(lhs, 15), mp.nstr(rhs, 15), lhs < rhs))
# A: a=1/2, beta=2: M = e^z (1 - 4z/3 + 4z^2/3), zero at 1/2 + i sqrt2/2
zA = mp.mpf(1)/2 + 1j*mp.sqrt(2)/2
report("A", mp.mpf(1)/2, 2, zA)
for z in [mp.mpc(0.3, 0.2), mp.mpc(-2, 1.5), mp.mpc(4, -3)]:
    q = mp.e**z*(1 - mp.mpf(4)/3*z + mp.mpf(4)/3*z**2)
    print("  A identity residual at", z, mp.nstr(abs(M(mp.mpf(1)/2, -mp.mpf(3)/2, z) - q), 5))
# C: a=3/2, beta=3: M = e^z (1 - 2z + 4z^2 + 8z^3/3)
roots = mp.polyroots([mp.mpf(8)/3, 4, -2, 1], maxsteps=200, extraprec=200)
for r in roots:
    if abs(mp.im(r)) > 1e-20:
        report("C", mp.mpf(3)/2, 3, r)
        p = lambda x: 1 - 2*x + 4*x**2 + mp.mpf(8)/3*x**3
        print("  C p(Re z)=", mp.nstr(p(mp.re(r)), 15))
        break
# B: a=2, beta=3/2 (b=1/2>0), numerical
report("B", 2, mp.mpf(3)/2, mp.mpc(-10, 5))
