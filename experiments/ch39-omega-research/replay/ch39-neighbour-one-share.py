"""允许各邻点重选共同标尺及第 0 条线分拆份额后的全域引力校准。"""
from flint import arb


def check(code, atoms, edge, rows, heat, moments, gauge, coefficients, initial, beta, logscale):
    pi = arb.pi()
    a, b, weight = edge
    common = atoms[2:]
    ends = [(a, weight), (b, weight)]
    panels = [[*common, end] for end in ends]
    C = arb(1)/500
    lo, hi = arb(".0068"), arb(".00702")
    span = lo.union(hi)

    def integral(c, panel, with_gauge=False):
        total = [arb(0) for _ in range(5 if with_gauge else 2)]
        for v, w in panel:
            x = c*v
            values = (gauge(rows, moments, heat, coefficients, x)[1:]
                      if with_gauge else moments(rows, heat, x))
            for j, val in enumerate(values):
                total[j] += w*val
        return total

    kprime, wprime = [], []
    for panel in panels:
        hterms = [(w*v*alpha/(2*pi), z/v) for v,w in panel for _,alpha,z in rows]
        slope = sum(h for h,_ in hterms)
        ceiling = C*(slope.max(0)+sum((-h).max(0)*(1-heat(2,y/C)) for h,y in hterms))
        floor = slope-sum(abs(h)*(1-(-y/C).exp()) for h,y in hterms)
        assert ceiling < arb(".353") and floor > arb("133.4")
        assert integral(lo,panel)[0] < 1 < integral(hi,panel)[0]
        dk = sum(h*(-y/span).exp() for h,y in hterms)
        dw = -sum(w*v*v*g*(-z/(v*span)).exp() for v,w in panel for g,_,z in rows)/(4*pi)
        kprime.append(dk)
        wprime.append(span*dw)

    def direction(c, with_gauge=False):
        left = integral(c,ends[:1],with_gauge)
        right = integral(c,ends[1:],with_gauge)
        return [l-r for l,r in zip(left,right,strict=True)]

    dk, dw = direction(span)
    assert dk > arb(".026")
    kc = kprime[0].union(kprime[1])
    assert kc > 0
    along = dw-wprime[0].union(wprime[1])*dk/kc
    assert along < arb("-7e-6")

    def at_newton(c, with_gauge=False):
        base = integral(c,panels[1],with_gauge)
        diff = direction(c,with_gauge)
        f = (1-base[0])/diff[0]
        return f, [u+f*d for u,d in zip(base,diff,strict=True)]

    brackets = {2345665:(".00699520267",".00699520270"),2345667:(".00699679978",".00699679981")}
    lower, upper = map(arb,brackets[code])
    assert at_newton(lower)[1][1] < arb("1e-122") < at_newton(upper)[1][1]
    f, values = at_newton(lower.union(upper),True)
    assert 0 < f < 1
    u = [v+a0-be*logscale for v,a0,be in zip(values[2:],initial,beta,strict=True)]
    g1, g2 = u[0]-u[2], u[1]-u[2]
    inv = g1-(beta[0]-beta[2])/(beta[1]-beta[2])*g2
    print("调标尺及一份额后的唯一引力解：零点、左端份额、两条规范差、消去能标的残差",
          code, *[v.str(12) for v in (f,g1,g2,inv)])
    assert inv < -19 and g1 > 220 and g2 > 83
