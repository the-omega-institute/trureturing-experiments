"""在规范均值与四项味窗下求端口 6 的无 11431 混合的最小读数偏差。"""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_mat, arb_poly, ctx, fmpq, fmpq_poly

root = Path(__file__).parent
single = run_path(str(root/"ch39-zero-single-port.py"))
sharp = run_path(str(root/"ch39-sharp-gauge-flavor-variance.py"))


def port_cost(include_tangent=False):
    q, *nums = single["polynomials"]()
    num = nums[6]
    off = sharp["study"]["gauge"]["OFFSET"]
    offset = fmpq(off.numerator, off.denominator)
    k = offset-fmpq(1, 4)
    den = fmpq_poly([1, -2*k])
    slope = num.derivative()*q-num*q.derivative()
    curvature = den*(slope.derivative()*q-2*slope*q.derivative())-4*k*slope*q
    q, p, slope, curvature = [arb_poly(poly) for poly in (q, num, slope, curvature)]
    domain = arb(-1).union(0)/200
    assert q(domain) > 0 and curvature(domain) < 0 and 0 not in slope(domain)
    assert acb.zeta_zero(1).imag > 10
    scale = acb.zeta_zero(11431).imag**2+arb(offset)
    def data(x):
        g = (1+x)/scale
        s = -g/(2*(1-arb(k)*g))
        factor = 1-2*arb(k)*s
        value = p(s)/q(s)
        first = -factor**2*slope(s)/(2*scale*q(s)**2)
        return value, first, s, factor
    center, tangent, s0, d0 = data(0)
    strength = -d0**3*curvature(s0)/(8*scale**2*q(s0)**3)
    def cost(x):
        value, first, _, _ = data(x)
        return (center+tangent*x-value)/strength, (tangent-first)/strength
    return (cost, strength, tangent) if include_tangent else (cost, strength)


def interval_distribution(rows, shifts, theta):
    # 在参数区间上同时包围取等权重和相位方程的唯一根。
    from flint import arb_poly
    coords, observations = zip(*rows)
    a, b, u, v, j = zip(*observations)
    S, lo, hi, tau = (sharp[name] for name in ("S","lo","hi","tau"))
    matrix = arb_mat([[1]*4, list(coords), [a0-S*(1-u0) for a0,u0 in zip(a,u)], list(u)])
    base = matrix.solve(arb_mat([[1],[0],[0],[lo]]))
    speed = matrix.solve(arb_mat([[0],[0],[0],[hi-lo]]))
    direction = matrix.solve(arb_mat([[0],[1],[0],[0]]))
    def dot(column):
        return sum(s*column[i,0] for i,s in enumerate(shifts))
    factor = theta/(1+theta*dot(direction))
    base = base-direction*(factor*dot(base))
    speed = speed-direction*(factor*dot(speed))
    weights = [arb_poly([base[i,0],speed[i,0]]) for i in range(4)]
    U = arb_poly([lo,hi-lo])
    B,V,J = [sum(w*y for w,y in zip(weights,col)) for col in (b,v,j)]
    P = (1-S)*U*(1-U)-(1-U)*V+(S-(1-S)*U)*B
    phase = J-(tau/2)*P
    box = arb(0).union(1)
    assert phase(0)*phase(1)<0 and 0 not in phase.derivative()(box)
    for _ in range(5):
        mid = box.mid()
        box = box.intersection(mid-phase(mid)/phase.derivative()(box))
    probs = [w(box) for w in weights]
    assert all(w>0 for w in probs) and J(box)>0
    assert arb(".041775")**2 < B(box)/(1-U(box)) < arb(".041885")**2
    return probs, [p(box) for p in (U,B,V,J)]


def certify_direction(theta, rows, costs, tails, scale):
    shifts = [scale*y for y in costs]
    probs, ref = interval_distribution(rows, shifts, theta)
    minimum = sum(w*y for w,y in zip(probs,costs))
    endpoint_duals = []
    for u in (sharp["lo"], sharp["hi"]):
        tests = sharp["phase_tests"](u, ref)
        matrix = arb_mat([[1, x+theta*shift, *tests(obs)]
                          for (x,obs),shift in zip(rows,shifts)])
        c = matrix.solve(arb_mat([[y] for y in costs]))
        assert matrix.det()>0 and c[2,0]>0 and c[3,0]>0
        assert (1-scale*theta*c[1,0])>0 and (1+scale*(1-theta)*c[1,0])>0
        assert c[0,0].overlaps(minimum)
        def linear(obs, c=c, tests=tests):
            a, b = tests(obs)
            return c[2,0]*a+c[3,0]*b
        def normalized(obs):
            a, b = tests(obs)
            return (c[2,0]/c[3,0])*a+b
        upper = c[3,0]*sharp["compare_torus"](arb(1), normalized)
        for n, x, value, dx in tails:
            speed = 1+scale*theta*dx
            assert speed>0
            assert (dx/speed > c[1,0] if n==8999 else dx/speed < c[1,0])
            assert value-c[1,0]*(x+theta*scale*value) > minimum+upper
        endpoint_duals.append((c[1,0],linear))
    return minimum, endpoint_duals


def verify():
    ctx.prec = 224
    cost, strength, tangent = port_cost(include_tangent=True)
    scale = -strength/tangent
    assert scale > 0
    rows = [sharp["row"](acb.zeta_zero(n)) for n in sharp["support"]]
    values = [cost(x)[0] for x, _ in rows]
    tails = []
    for n in (8999,15001):
        x, _ = sharp["row"](acb.zeta_zero(n))
        tails.append((n, x, *cost(x)))
    for theta in (arb(0),arb(1)/2,arb(1)):
        minimum, _ = certify_direction(theta,rows,values,tails,scale)
        print("Budget share, relative threshold, port deficit:",
              theta.str(3), (scale*minimum).str(20),
              (-tangent*(1-theta)*scale*minimum).str(20), flush=True)
    # 环面仿射界 <0.06 由 ch39-gauge-flavor-variance.py 的 torus_bound 给出。
    intervals = [arb(k).union(k+1)/128 for k in range(128)]
    checks = [(theta, *certify_direction(theta,rows,values,tails,scale)) for theta in intervals]
    floor = arb(1)
    for start in range(9000,15001,256):
        count = min(256,15001-start)
        zeros = acb.zeta_zeros(start,count)
        assert len(zeros)==count
        for n, zero in enumerate(zeros,start):
            if n==11431 or n in sharp["support"]:
                continue
            x, obs = sharp["row"](zero)
            value, _ = cost(x)
            for theta, minimum, duals in checks:
                coordinate = x+theta*scale*value
                for slope, linear in duals:
                    margin = value-minimum-slope*coordinate-linear(obs)
                    assert margin > 0, (n, str(theta), str(margin))
                    floor = min(floor,margin.lower())
    print("Uniform strict slack over all directions:",floor.str(12),flush=True)


if __name__ == "__main__":
    verify()
