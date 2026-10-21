import sympy as sp
s,e=sp.symbols('s epsilon',positive=True)
def coeffs(g,n=3,order=6):
    ser=sp.series(g*s**n,s,0,3).removeO()
    c=[sp.simplify(ser.coeff(s,k)) for k in range(3)]
    A0=c[0]/sp.factorial(n-1); A1=c[1]/sp.factorial(n-2)
    return sp.simplify(A0),sp.simplify(A1)
def futaki(gzeta,geta,n=3):
    A0,A1=coeffs(gzeta,n); A0e,A1e=coeffs(geta,n)
    B0=-sp.diff(A0e,e).subs(e,0)/n; B1=-sp.diff(A1e,e).subs(e,0)/(n-1)
    return sp.simplify(A1/A0*B0-B1), A0, A1
t=sp.exp(-s)
cases={
 'C3 orbifold check: C^3 (R-charges 2/3), weight x':((1/(1-t**sp.Rational(2,3))**3), 1/((1-t**(sp.Rational(2,3)+e))*(1-t**sp.Rational(2,3))**2)),
 'conifold, generator R=1':((1-t**2)/(1-t)**4, (1-t**2)/((1-t**(1+e))*(1-t)**3)),
}
for r in [2,3]:
    # xy = z^r w^r with R(z)=R(w)=1, R(x)=R(y)=r: g = (1-t^{2r})/((1-t^r)^2 (1-t)^2)
    cases[f'xy=z^{r}w^{r}, generator z (R=1)']=((1-t**(2*r))/((1-t**r)**2*(1-t)**2),(1-t**(2*r))/((1-t**r)**2*(1-t**(1+e))*(1-t)))
for name,(g0,ge) in cases.items():
    F,A0,A1=futaki(g0,ge)
    V=sp.Rational(16,27)*A0
    print(f'{name}: A0={A0} A1={A1} F={F} (27/64)V={sp.Rational(27,64)*V} F>=bound? {sp.simplify(F-sp.Rational(27,64)*V)>=0}')
