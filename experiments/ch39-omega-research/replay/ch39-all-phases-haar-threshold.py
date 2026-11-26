"""正谱矩与四维 Haar 求积给任意单位边相位的一圈匹配约束。"""
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import acb, acb_mat, arb, arb_poly, ctx


def flat_walk(shell, reduction):
    matrix = acb_mat(19, 19)
    for i in range(19):
        matrix[i,i] = 1
    pole = shell["soft"]["pole"]
    tuned = shell["soft"]["star"]["orbit"]
    for i, j, _, _ in shell["geo"]["edges"]:
        spoke = pole in (i,j) and (i in tuned or j in tuned)
        step = (arb(1)-reduction if spoke else arb(1))/36
        matrix[i,i] -= step
        matrix[j,j] -= step
        matrix[i,j] += step
        matrix[j,i] += step
    return matrix


def envelopes(steps, lower, upper, period):
    moment = lambda matrix: sum(matrix[i,i].real for i in range(19))
    tau = moment(lower**steps)
    a = (moment(lower**2)/moment(upper))**steps
    raw = list((upper**steps).charpoly())
    poly = arb_poly([v.real for v in reversed(raw)])
    derivative = poly.derivative()
    arrays = [np.empty(period) for _ in range(4)]
    for i in range(period):
        if i in (0,period//2):
            u, du = arb(int(i != 0)), arb(0)
        else:
            sine, cosine = (arb.pi()*i/period).sin_cos()
            u, du = sine*sine, sine*cosine
        at_pole = i == period//2
        p = arb(0) if at_pole else poly(u)
        low_slope = tau*du/(1-a*u)
        high_slope = arb(0) if at_pole else -du*derivative(u)/p
        mid = (high_slope+low_slope)/2
        width = abs(high_slope-low_slope)/2
        top = (1-a*u)**(6*tau/a)
        values = (p**6 if not at_pole else arb(0), top, mid, width)
        for j, (arr, x) in enumerate(zip(arrays, values)):
            arr[i] = float(x)
            tolerance = abs(x)*arb("1e-15")
            if j >= 2:
                tolerance += arb("1e-45")
            assert x == 0 or abs(arb(arr[i])-x) < tolerance
    assert max(np.max(np.abs(x)) for x in arrays[2:]) < 1100
    return arrays


def verify():
    directory = Path(__file__).parent
    shell = run_path(str(directory/"ch39-three-ray-chiral-boundary.py"))
    grid = run_path(str(directory/"ch39-three-mode-haar-grid.py"))
    full = run_path(str(directory/"ch39-full-shell-haar-grid.py"))
    ctx.prec = 320
    max_reduction = arb("1e-8")
    offset = shell["soft"]["point"]["study"]["offset"]
    assert arb(23)/25*arb(2)**(-20)/(acb.zeta_zero(1).imag**2+offset) < max_reduction
    lower, upper = flat_walk(shell,0), flat_walk(shell,max_reduction)
    assert sum(upper[i,i].real for i in range(19)) < arb("11.501")
    _,_,weyl = grid["tables"](840)
    alias = full["alias_bound"](arb("11.501"))
    for steps in (1,2):
        dlow, dhigh, mid, width = envelopes(steps,lower,upper,840)
        bounds = []
        for floor, density in ((True,dlow),(False,dhigh)):
            def rms(slope, absolute):
                factors = dict.fromkeys((1,2,3,4,6),(density,slope))
                z, norms = grid["grid_moments"](840,210,210,factors,weyl,absolute)
                return [z.sqrt()*n for n in norms]
            center, errors = rms(mid,False), rms(width,True)
            if floor:
                assert all(c>e for c,e in zip(center,errors))
            sign = -1 if floor else 1
            bounds.append([(c+sign*e)**2 + sign*alias*k*k*n for c,e,k,n in
                           zip(center,errors,(684,152,152),(1,1,3))])
        y,w,c = [lo.lower().union(hi.upper()) for lo,hi in zip(*bounds)]
        # Eliminate the unknown UV scale and logarithmic running time with
        # the normal of the plane through the low-energy data and the beta vector.
        a1 = 3*(1-arb(".23129"))*arb("127.930")/5
        a2, a3 = arb(".23129")*arb("127.930"), 1/arb(".1180")
        b1,b2,b3 = arb(41)/10,-arb(19)/6,-7
        n_y = (a2*b3-a3*b2)/(a1*b2-a2*b1)
        n_w = (a3*b1-a1*b3)/(a1*b2-a2*b1)
        hy, hc = y/(5*w), c/(4*w)
        gap = n_y*hy+n_w+hc
        # With the electroweak inputs fixed, alpha_s decreases separately
        # with both ratios here. Bounding the four corners also works by
        # separate linear-fractional monotonicity and positive denominators.
        strong = [grid["one_loop"]([60*r,s])[1]
                  for r in (hy.lower(),hy.upper()) for s in (hc.lower(),hc.upper())]
        assert all(a > 0 for a in strong)
        if steps == 1:
            assert arb(".1322") < gap < arb(".1324")
            assert all(arb(".07100") < a < arb(".07105") for a in strong)
        else:
            assert arb(".1157") < gap < arb(".1517")
            assert all(arb(".06730") < a < arb(".07464") for a in strong)
        print("步数、超荷/六十倍弱、色/弱、一圈投影缺口",
              steps, hy.lower().str(9),hy.upper().str(9),
              hc.lower().str(9),hc.upper().str(9),
              gap.lower().str(9),gap.upper().str(9),flush=True)


if __name__ == "__main__":
    verify()
