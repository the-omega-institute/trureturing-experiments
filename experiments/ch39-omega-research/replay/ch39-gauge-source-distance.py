"""Bracket the smallest maximum relative displacement of the low gauge inputs."""
from pathlib import Path
from runpy import run_path

from flint import arb

floating = run_path(str(Path(__file__).with_name("ch39-floating-clock-ground.py")))
gauge = floating["gauge"]
shortest = gauge["shortest"]
model = shortest["model"]
volume = gauge["events"]["mixing"]["volume"]
q = shortest["q"]


def verify():
    offsets = [q(x) for x in (".000594660007738019","-.000647357405682778",".000647357405682778")]
    factors = [1+h for h in offsets]
    ref_g = [c*f for c,f in zip(model["centers"],factors)]
    ref_y = model["yukawa_centers"].copy()
    ref_y[2] *= q(".99")
    time = q("26.816785428385903")
    state = [u*u for u in ref_g+ref_y]
    for _ in range(128):
        state = shortest["reference"]["advance"](state,time/128,14)
    common = model["target"](11431,time)
    residuals = [1/(r*u)-common for r,u in zip(shortest["ratios"],state)]
    mix = [a+b for a,b in gauge["events"]["mixing"]["displacement"](27)]
    errors = [abs(e)+d for e,d in zip(residuals,mix)]
    assert max(e.upper() for e in errors) < q("1.01e-7") < q("3.18e-7")
    # Move the two limiting inputs inward and balance the resulting residuals
    # by a small common scale shift. The first source factor stays fixed.
    ppm = q("647.291")/10**6
    close_g = [ref_g[0],model["centers"][1]*(1-ppm),model["centers"][2]*(1+ppm)]
    close_time = time+q(".0000000140")
    state = [u*u for u in close_g+ref_y]
    for _ in range(128):
        state = shortest["reference"]["advance"](state,close_time/128,14)
    close_errors = [abs(1/(r*u)-model["target"](11431,close_time))+d
                    for r,u,d in zip(shortest["ratios"],state[:3],mix,strict=True)]
    assert max(e.upper() for e in close_errors) < q("3.18e-7")
    print("647.291 ppm all-mixing trial errors",*(e.upper().str(8) for e in close_errors))
    xl,xh,logs = volume["input_bounds"]()
    jac,_ = volume["variations"]([l.union(h) for l,h in zip(xl,xh)],logs,refinements=2)
    speeds = floating["flow_bounds"]()
    anchor = q("429/16")
    low, high = [model["target"](11431,anchor+s*q("1/8")) for s in (1,-1)]
    target_range = low.union(high)
    beta = 1+(speeds[2]-speeds[1])/(2*target_range+speeds[1])
    assert beta > 1
    row = [jac[2,j]/shortest["ratios"][2]-beta*jac[1,j]/shortest["ratios"][1]
           for j in range(12)]
    assert row[1]<0<row[2]
    y_changes = [y-2*ref.log() for y,ref in zip(logs,ref_y)]
    # Ordinary wide ball arithmetic can overestimate products on one side;
    # use the four endpoint products for a sign-aware lower bound.
    def sharp_floor(ppm):
        amount = q(ppm)/10**6
        g_changes = [(1/(c*(1+amount))**2-1/ref**2,1/(c*(1-amount))**2-1/ref**2)
                     for c,ref in zip(model["centers"],ref_g)]
        changes = g_changes+[(x.lower(),x.upper()) for x in y_changes]
        v = (residuals[2]-beta*residuals[1]).lower()
        for a,(lo,hi) in zip(row,changes):
            v += min((endpoint*x).lower() for endpoint in (a.lower(),a.upper()) for x in (lo,hi))
        return (v-mix[2]-beta*mix[1]).lower()
    tol = arb(2)**-21+arb("1.3e-16")
    floor = sharp_floor("647.17")
    assert floor > tol*(1+beta)
    print("excluded radius 647.17 ppm, separation",floor.str(9))
    print("explicit witness",max(abs(v).upper() for v in offsets).str(15),
          "all mixing errors below",max(v.upper() for v in errors).str(9))


if __name__ == "__main__":
    verify()
