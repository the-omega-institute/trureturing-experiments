"""Uniform error of the rational gauge and target registers over all finite input codes."""
from math import factorial
from pathlib import Path
from runpy import run_path

from flint import acb, arb, ctx

root = Path(__file__).parent
reference = run_path(str(root/"ch39-grid-gauge-readout.py"))
volume, model = reference["volume"], reference["model"]
exact, q, running = reference["rational"], reference["q"], reference["running"]


def flow_error():
    low,high,logs = volume["input_bounds"]()
    domain = []
    def observe(step,u,y,tr):
        if domain:
            domain[:] = [a.union(b) for a,b in zip(domain,u+y)]
        else:
            domain.extend(u+y)
    model["evolve"]([a.union(b) for a,b in zip(low,high)],logs,q(33),8192,
                    diagonal_quarks=True,observe_flavors=observe)
    padded = [v+arb(0,"1e-7") for v in domain]
    u,y = padded[:3],padded[3:]
    tr = [sum(y[i:i+3]) for i in (0,3,6)]
    grid = arb(1)/exact["SCALE"]
    k = q(exact["INVERSE_LOOP"])
    assert 0 < 1/model["loop"]-k < grid
    k_range = k.union(1/model["loop"])
    # The following row sums and parameter derivatives hold on the padded box.
    state_rates, parameter_rates = [], []
    for i in range(3):
        b,d,c = [running[key][i] for key in ("B","D","C")]
        corr = sum(a*x-z*t for a,z,x,t in zip(b,d,u,tr))
        beta = running["BETA"][i]
        state_rates.append(4*abs(u[i]*(beta+corr*k_range)*k_range)
            +2*u[i]**2*k_range**2*(sum(map(abs,b))+3*sum(map(abs,d))))
        parameter_rates.append(abs(2*u[i]**2*(beta+2*k_range*corr)))
        gy = sum(a*x for a,x in zip(c,u))
        for j in range(3):
            a = y[3*i+j]
            mix = y[3*(1-i)+j] if i<2 else arb(0)
            field = 3*(a-mix)+2*(3*(tr[0]+tr[1])+tr[2]-gy)
            state_rates.append(k_range*(abs(field)+abs(a)*(48+2*sum(map(abs,c)))))
            parameter_rates.append(abs(a*field))
    lipschitz = max(v.upper() for v in state_rates)
    perturbation = grid*max(v.upper() for v in parameter_rates)
    assert lipschitz < q(".457")
    terms = exact["taylor"](padded,exact["DEGREE"]+1,running,k)
    remainder = max(abs(v).upper() for v in terms[-1])*(q(33)/exact["STEPS"])**(exact["DEGREE"]+1)
    bound = (33*lipschitz).exp()*(grid+33*perturbation+exact["STEPS"]*(remainder+grid))
    assert bound < q("1e-15")
    # The same estimate bounds the within-step separation before each local
    # truncation; by first exit, the rational-parameter comparison flow stays in the pad.
    assert min(v.lower() for v in u) > q(".2")
    return bound/(min(model["ratios"])*q(".2")**2)


def verify():
    ctx.prec = 176
    for key in ("BETA","C","B","D"):
        rows = [running[key]] if key == "BETA" else running[key]
        data = exact["COEFFICIENTS"][key]
        expected = [data] if key == "BETA" else data
        assert all(x.contains(q(a).mid()) or (x-q(a)).contains(0)
                   for row,items in zip(rows,expected,strict=True)
                   for x,a in zip(row,items,strict=True))
    assert all((v-q(t)).contains(0) for v,t in zip(model["centers"]+model["ratios"]+model["yukawa_centers"],
               exact["GAUGE"]+exact["RATIOS"]+exact["YUKAWA"],strict=True))
    gauge_error = flow_error()
    fixed = model["fixed"]
    common = (fixed["geo"]["BOOK_PLANCK"]/1000)**2*3*arb(5).sqrt()/(4*arb.pi())*arb(-46).exp()
    register = exact["register"]()
    assert list(register) == list(range(10815,12093))
    for n,z in enumerate(acb.zeta_zeros(10815,1278),10815):
        value = common/(z.imag**2+q(fixed["OFFSET"]))
        assert (value*2**96).floor().unique_fmpz() == register[n]
        assert 0 < value < 10000
    # Even truncated alternating Taylor polynomial on [-20/128,0] lies in [0,1].
    # Squaring on [0,1] has Lipschitz constant 2, each truncation loses < 2^-128.
    clock_error = 128*(q(20)/128)**21/factorial(21)+q(127)/2**128
    target_error = q(1)/2**96 + 10000*clock_error + q(1)/2**128
    total = gauge_error+target_error+q(1)/2**43
    assert total < q(1)/2**40
    print("All-code absolute gauge, target and rounded residual error bounds",
          gauge_error.upper().str(8),target_error.upper().str(8),total.upper().str(8))


if __name__ == "__main__":
    verify()
