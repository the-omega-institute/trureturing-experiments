"""Novaes 2023 (arXiv:2305.18235, J. Stat. Phys. 190, 182), conjecture after Eq. (53): Eqs. (54) and (55).

Model: the beta = 2 law of the repository's StrongAbsorption.Foundations.law, density proportional to
exp(-M(sum x + sum y)) Delta(x)^2 Delta(y)^2 prod_{i,j} (gamma + x_i - y_j) on x_i >= 0, 0 <= y_j <= gamma, with proper
delay times q_i = 1/(gamma + x_i). With t = 1/gamma, prod (gamma + x_i - y_j) = gamma^(M^2) prod (1 + t(x_i - y_j)); the
cutoff y_j <= gamma changes polynomial expectations only by exp(-M gamma / 2) terms, so the expansion in t is
E0[f D_t] / E0[D_t] with E0 the product of two Laguerre (beta = 2) ensembles. Exact rational arithmetic (sympy).
Checks: the t-series of gamma^n (1/M)<Tr Q^n> through t^3 equals Eq. (54) for 1 <= M <= 3, 1 <= n <= 4 (n <= 3 at M = 3),
and the t-series of gamma^n <(Tr Q / M)^n> through t^2 equals Eq. (55) for 1 <= M <= 3, 1 <= n <= 3; Eq. (54) at n = 2
reproduces the source's Eq. (52).
"""
import sympy as sp, sys
from math import factorial
t=sp.symbols('t')
def E(M, f, order):
    xs=sp.symbols(f'x0:{M}'); ys=sp.symbols(f'y0:{M}')
    vand=lambda v: sp.prod([(v[i]-v[j])**2 for i in range(M) for j in range(i+1,M)])
    D=sp.prod([1+t*(xi-yj) for xi in xs for yj in ys]); base=vand(xs)*vand(ys)
    def E0(poly):
        P=sp.Poly(sp.expand(poly*base), *xs, *ys); tot=0
        for mon,coef in P.terms():
            val=1
            for a in mon: val*=sp.Rational(factorial(a), M**(a+1))
            tot+=coef*val
        return sp.expand(tot)
    ff=f(xs)
    num=E0(sp.expand(sp.series(ff*D,t,0,order).removeO())); den=E0(D)
    return sp.expand(sp.series(num/den,t,0,order).removeO())
ok=True
for M in (1,2,3):
    for n in (1,2,3):
        r54=E(M, lambda xs: sum(sp.series((1+t*x)**(-n),t,0,4).removeO() for x in xs)/M, 4)
        c54=1-n*t+n**2*t**2-sp.Rational((n+1)*(n*(5*n-2)*M**2+n*(n+2)),6*M**2)*t**3
        r55=E(M, lambda xs: sp.expand(sp.series((sum(sp.series((1+t*x)**(-1),t,0,3).removeO() for x in xs)/M)**n,t,0,3).removeO()), 3)
        c55=1-n*t+sp.Rational(n*((n+1)*M**2+(n-1)),2*M**2)*t**2
        d54=sp.expand(r54-c54); d55=sp.expand(r55-c55); ok = ok and d54==0 and d55==0
        print(f"M={M} n={n} (54)series={r54} diff={d54}; (55)series={r55} diff={d55}", flush=True)
    if M < 3:
        n=4
        r54=E(M, lambda xs: sum(sp.series((1+t*x)**(-n),t,0,4).removeO() for x in xs)/M, 4)
        c54=1-n*t+n**2*t**2-sp.Rational((n+1)*(n*(5*n-2)*M**2+n*(n+2)),6*M**2)*t**3
        d54=sp.expand(r54-c54); ok = ok and d54==0
        print(f"M={M} n={n} (54)series={r54} diff={d54}", flush=True)
eq52 = all(sp.expand(1-2*t+4*t**2-sp.Rational((2+1)*(2*(5*2-2)*M**2+2*(2+2)),6*M**2)*t**3 - (1-2*t+4*t**2-sp.Rational(4*(2*M**2+1),M**2)*t**3))==0 for M in range(1,20))
ok = ok and eq52
print("Eq. (54) at n = 2 equals Eq. (52) for M < 20:", eq52)
print("RESULT", "PASS" if ok else "FAIL")
raise SystemExit(0 if ok else 1)
