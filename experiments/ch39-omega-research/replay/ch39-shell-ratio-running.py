"""Transverse running of the family shell ratios at fixed one-loop matching."""
from pathlib import Path
from runpy import run_path


from flint import arb, arb_series, ctx
from sympy import cancel, symbols

TEN_ONE, TEN_SELF, TEN_CROSS = run_path(str(Path(__file__).with_name("ch39-two-loop-product-beta.py")))["target_inverse_coefficients"]()


def check_identities():
    e, s = symbols("e s", positive=True)
    x = [e, 1, s]
    y = [e*t for t in x]
    dx = [10*t+sum(x)+5*(y[a]+y[(a+1)%3]) for a,t in enumerate(x)]
    dy = [10*t+sum(y)+5*(x[a]+x[a-1]) for a,t in enumerate(y)]
    slopes = [b-a for a,b in zip(dx,dy)]
    average = sum(slopes)/3
    variance = sum((val-average)**2 for val in slopes)/3
    assert cancel(average-(e-1)*sum(x)) == 0
    ratio = variance/average**2
    floor = 25*(1-e+e**2)/(6*(1+e+e**2))
    remainder = (25*(1-e+e**2)*(s*(1+e)-1-e**2)**2
                 /(2*(1-e)**2*(1+e+s)**2*(1+e+e**2)))
    assert cancel(ratio-floor-remainder) == 0


def verify_flow():
    ext = run_path(str(Path(__file__).with_name("ch39-selected-extended-flow.py")))
    a = ext["a"]
    v, (g,h,old_x,old_y,f,j) = ext["initial"]()
    e = old_y/old_x
    floor = (25*(1-e+e**2)/(6*(1+e+e**2))).sqrt()
    assert a("2.0412") < floor < a("2.0413")
    x0 = [9*e,arb(9),arb(9)]
    state = [g,h,*x0,*[e*t for t in x0],f,j]
    pi = arb.pi()
    u = 1/(8*pi*pi)

    def rhs(state):
        g,h = state[:2]
        x,y = state[2:5],state[5:8]
        f,j = state[-2:]
        z,w = 1/g,1/h
        sx,sy = sum(x),sum(y)
        guard = a("54/5")*z+a("297/20")*w
        dg = u*(a("2/3")-u/2*(a("8192/15")*z+594*w-15*(sx+sy)-a("12/5")*f-a("54/5")*j))
        dh = u*(TEN_ONE+u/2*(TEN_SELF*w-TEN_CROSS*z+5*(sx+sy)))
        dx = [u*t*(10*t+sx+5*(y[i]+y[(i+1)%3])-guard) for i,t in enumerate(x)]
        dy = [u*t*(10*t+sy+5*(x[i]+x[i-1])-guard) for i,t in enumerate(y)]
        return [dg,dh,*dx,*dy,u*f*(a("16/5")*f+3*j-a("72/5")*z),
                u*j*(f+a("44/5")*j-a("108/5")*z)]

    def taylor(box, order):
        cols = [[value] for value in box]
        for i in range(order):
            slope = rhs([arb_series(col,prec=i+1) for col in cols])
            for col, deriv in zip(cols,slope):
                col.append(deriv[i]/(i+1))
        return cols

    step = a("1/1024")
    speeds = [1,1,a("1e-7"),32,32,a("1e-15"),a("1e-7"),a("1e-7"),32,32]
    margins = [(-s*step).union(s*step) for s in speeds]
    for i in range(800):
        box = [t+m for t,m in zip(state,margins)]
        deriv = rhs(box)
        assert all(abs(d)<s for d,s in zip(deriv,speeds))
        assert deriv[3]>0 and deriv[4]>0
        assert all(0 < box[j] < 4*pi for j in (2,5,6,7,8,9))
        assert 1/(4*pi*box[0])<1/4 and 1/(4*pi*box[1])<1/4
        coeff, rem = taylor(state,3), taylor(box,4)
        nxt = []
        for col,last in zip(coeff,rem):
            value = last[-1]
            for c in reversed(col):
                value = c+step*value
            nxt.append(value)
        if max(nxt[3].lower(),nxt[4].lower())>4*pi:
            assert state[3] < 4*pi and state[4] < 4*pi
            scale = v*((i*step).exp().union(((i+1)*step).exp()))
            assert a("8.35e14") < scale < a("8.37e14")
            print("Log step", i, "to", i+1, "over 1024")
            print("Family ratios relative to the shell at the last step start",
                  [(state[5+j]/(state[2+j]*e)).sqrt().str(8) for j in range(3)])
            print("First bridge squared Yukawa at 4*pi / GeV",scale.str(8))
            return scale
        assert nxt[3]<4*pi and nxt[4]<4*pi
        state = nxt
    raise AssertionError("Missing crossing")


if __name__ == "__main__":
    ctx.prec = 150
    check_identities()
    verify_flow()
