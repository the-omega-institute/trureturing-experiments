#!/usr/bin/env python3
"""Exact algebra checks for Chebotarev arXiv:2611.02249v1, Section 5.

No floating point is used. These are preregistration checks, not a Lean proof.
"""
import sympy as sp

u, y, a, b, al, be = sp.symbols("u y a b al be")
w = sp.symbols("w", real=True)


def zero(expr):
    assert sp.cancel(expr) == 0, expr


def universal_identities():
    z = sp.symbols("z")
    block = sp.Matrix([[a + w*al/2, w*be/2], [w*al/2, b + w*be/2]])
    characteristic = (z-a)*(z-b) - w*(be*(z-a)+al*(z-b))/2
    zero(block.charpoly(z).as_expr() - characteristic)
    discriminant = (b-a)**2 + w*(b-a)*(be-al) + w**2*(al+be)**2/4
    zero(sp.discriminant(characteristic, z)-discriminant)
    h, c, s, S, C = sp.symbols("h c s S C", real=True)
    # c=cos(x), s=sin(x), S=sin(theta), C=cos(theta).
    F = y*y - (2*h*(c+sp.I*s)+w*(C+sp.I*S))*y + 2*sp.I*w*h*S*(c+sp.I*s)
    t = 2*w*h*c*S/(2*h*c+w*C)
    value = sp.expand_complex(F.subs(y, sp.I*t))
    zero(sp.im(value))
    expected = 2*h*w*w*S*C*(S*c-C*s)*(w-1)/(2*h*c+w*C)**2
    # The only trigonometric relation required is 2h sin(x)cos(x)=sin(theta-x).
    relation = 2*h*c*s-S*c+C*s
    zero((sp.re(value)-expected)*(2*h*c+w*C)**2 + 2*C*S*h*w*w*relation)
    zero((sp.re(value)-expected).subs(h, (S*c-C*s)/(2*c*s)))
    y1 = sp.I*S
    y2 = -sp.I*S*C*(c-sp.I*s)/(2*h)
    # The w^2 coefficient vanishes under c^2+s^2=1.
    coeff2 = sp.expand(F.subs(y, y1*w+y2*w*w)).coeff(w, 2)
    zero(coeff2 - sp.I*S*C*(c*c+s*s-1))
    zero(sp.re(y2) + S*C*s/(2*h))
    # Re R >= 1-v^2 >= 3/4 for 0<=v<=1/2 and cos((m+2)x)>=0.
    v, A, B = sp.symbols("v A B", real=True)
    real_R = 1+2*v*A+v*v*B
    zero(real_R-(1-v*v)-(2*v*A+v*v*(B+1)))
    zero(1-v*v-sp.Rational(3, 4)-(sp.Rational(1, 2)-v)*(sp.Rational(1, 2)+v))
    print("UNIVERSAL: quadratic, discriminant, crossing identity, departure coefficient, Re R bound OK")


def small_cases():
    cases = rows = 0
    for n in range(8, 41, 2):
        phi = sp.Poly(sp.cyclotomic_poly(n, u), u, domain=sp.QQ)

        def reduce(expr):
            num, den = sp.fraction(sp.cancel(expr))
            # Rational functions are zero iff their numerator vanishes modulo Phi_n;
            # also require that the denominator is invertible in Q(mu).
            pn = sp.Poly(num, u, domain=sp.QQ.frac_field(w))
            pd = sp.Poly(den, u, domain=sp.QQ.frac_field(w))
            assert sp.gcd(pd, phi).degree() == 0
            return pn.rem(phi).as_expr()

        def equal(expr):
            assert reduce(expr) == 0, (n, m, expr)

        for m in range(2, n, 2):
            if 2*(m+1) >= n:
                continue
            cases += 1
            p = m+1
            ka = sum(1-u**s for s in range(1, m+1))/n
            kb = sum(1-(-u)**s for s in range(1, m+1))/n
            alpha, beta = (1-u**p)/n, (1+u**p)/n
            q = 2*sum(1-u**(2*j) for j in range(1, m//2+1))/n
            equal(q-ka-u*(1-u**m)/(n*(1+u)))
            equal(ka+kb+sp.Rational(1,n)-q-sp.Rational(p,n))
            # Proposition 10: apply the original row formula to BOTH Fourier modes.
            for k in range(n):
                def fa(j):
                    return u**(j % n)
                def fb(j):
                    return (-1)**(j % n)*fa(j)
                for mode, rhs in [
                    (fa, (ka+w*alpha/2)*fa(k)+(w*alpha/2)*fb(k)),
                    (fb, (w*beta/2)*fa(k)+(kb+w*beta/2)*fb(k)),
                ]:
                    lhs = (sum(mode(k)-mode(k+s) for s in range(1,m+1))
                           +(w if k % 2 == 0 else 0)*(mode(k)-mode(k+p)))/n
                    equal(lhs-rhs)
                # Proposition 10 endpoint vector x_(2i)=x_(2i+1)=mu^(2i).
                def pair(j):
                    j %= n
                    return u**(2*(j//2))
                lhs = (sum(pair(k)-pair(k+s) for s in range(1,m+1))
                       +(1 if k % 2 == 0 else 0)*(pair(k)-pair(k+p)))/n
                equal(lhs-q*pair(k))
                rows += 1
            z = sp.symbols("z")
            P = (z-ka)*(z-kb)-w*(beta*(z-ka)+alpha*(z-kb))/2
            equal(P.subs({z:q,w:1}))
            equal(P.subs({z:sp.Rational(p,n),w:1}))
            # R defined by normalizing the actual discriminant by (kb-ka)^2.
            D = (kb-ka)**2+2*w*(kb-ka)*u**p/n+w*w/n**2
            hphase = u*(1-u**m)/(1-u**2)  # h exp(-imx)
            equal(kb-ka-2*hphase/n)
            equal(D-(2*hphase/n)**2*(1+w*u**p/hphase+w*w/(4*hphase**2)))
        print(f"CYCLOTOMIC: n={n}, exact restriction and endpoint checks OK")
    print(f"EXACT_SMALL_CASES: {cases} admissible (n,m), {rows} original-matrix rows, n<=40")


if __name__ == "__main__":
    print("sympy", sp.__version__, "; arithmetic: QQ and exact polynomial identities")
    universal_identities()
    small_cases()
    print("ALL_CHECKS_OK")
