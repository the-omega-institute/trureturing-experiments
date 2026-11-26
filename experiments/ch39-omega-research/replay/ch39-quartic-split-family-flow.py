"""Certified two-loop gauge and six one-loop Yukawa trajectories for split families."""
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path

from flint import arb, arb_series, ctx

ROOT = Path(__file__).parent
TEN_ONE, TEN_SELF, TEN_CROSS = run_path(str(ROOT/"ch39-two-loop-product-beta.py"))["target_inverse_coefficients"]()
PLANCK = run_path(str(ROOT/"ch39-quartic-planck-branches.py"))
RATIO = PLANCK["RATIO"]
a = RATIO["a"]


def initial_states():
    products = RATIO["product_polynomials"]()
    powers = PLANCK["mass_exponents"](products)
    scores = RATIO["thresholds"]
    match = RATIO["MODEL"]["matched"]
    k, large = a("159/100"), a("333/200")
    for n,height,eps in RATIO["high_zero_readings"]():
        low,v,_ = match(scores(products,arb(1)),eps,arb(1))
        anchor = v/low*(2*arb.pi()*(height**2+a("1/4")+30+a("2/15"))/(3*arb(5).sqrt()/2)).sqrt()
        def mass(r):
            return anchor*PLANCK["mass_ratio"](r,powers)
        bounds = PLANCK["enclose_root"](mass, PLANCK["TARGET"], F(".45"), F(".46"), False)
        ratios = [("center", a(bounds[0]).union(a(bounds[1])))]
        if n == 11485:
            rational = a("24/53")
            assert abs(mass(rational)-PLANCK["TARGET"]) < PLANCK["HALF_UNIT"]
            ratios.append(("rational", rational))
        for kind,r in ratios:
            low,v,lx = match(scores(products,r),eps,k,r)
            mean = lx.exp()
            small = mean**3/large**2
            y = eps*mean
            assert min(large,small,y)*v > 1000
            state = [1/(low**2*(1+3/k**2)),1/(low**2*(3+k**2)),
                     large**2,large**2,small**2,y*y,y*y,y*y]
            yield n, kind, v, state


def full_flow(state, duration):
    pi = arb.pi()
    unit = 1/(8*pi**2)
    b5, diag5, diag10 = a("10/3"),a("6796/15"),TEN_SELF
    cg5, cg10 = a("54/5"),a("297/20")
    def rhs(z):
        g,h = z[:2]
        x,y = z[2:5],z[5:]
        sx,sy = sum(x),sum(y)
        gauge = cg5/g+cg10/h
        dg = unit*(b5-(diag5/g+594/h-15*(sx+sy))*unit/2)
        dh = unit*(TEN_ONE+(diag10/h-TEN_CROSS/g+5*(sx+sy))*unit/2)
        dx = [unit*xi*(10*xi+sx+5*(y[i]+y[(i+1)%3])-gauge) for i,xi in enumerate(x)]
        dy = [unit*yi*(10*yi+sy+5*(x[i]+x[(i-1)%3])-gauge) for i,yi in enumerate(y)]
        return [dg,dh,*dx,*dy]

    def jet(box, degree):
        seqs = [[value] for value in box]
        for j in range(degree):
            derivatives = rhs([arb_series(s,prec=j+1) for s in seqs])
            for seq, slope in zip(seqs,derivatives):
                seq.append(slope[j]/(j+1))
        return seqs

    step = a("1/128")
    speeds = [2,2,32,32,32,2,2,2]
    margins = [(-s*step).union(s*step) for s in speeds]
    for i in range(int((duration/step).ceil().unique_fmpz())):
        tube = [s+m for s,m in zip(state,margins)]
        slope = rhs(tube)
        assert all(abs(ds)<s for ds,s in zip(slope,speeds))
        assert slope[0]<0<slope[1]
        # H is increasing, so its starting endpoint is the lower bound for the step.
        assert 1/(4*pi*state[1]) < a("1/4") and 1/(4*pi*tube[0]) < a("1/4")
        assert all(value < 4*pi for value in tube[2:])
        if (i+1)*step >= duration.upper():
            dt = duration-i*step
            end = [s+dt*f for s,f in zip(state,slope)]
            return 1/(4*pi*end[0]),1/(4*pi*end[1])
        exact, remainder = jet(state,4), jet(tube,5)
        state = []
        for coeff, last in zip(exact,remainder):
            value = last[-1]
            for c in reversed(coeff):
                value = c+step*value
            state.append(value)
    raise AssertionError("Target time was not enclosed")


def verify():
    ctx.prec = 160
    for n,kind,v,state in initial_states():
        alpha5,alpha10 = full_flow(state,(a(PLANCK["TARGET"])/v).log())
        assert a(".22") < alpha5 < a(".25") and a(".017") < alpha10 < a(".019")
        print(n,kind,"Planck alpha5,alpha10",alpha5.str(8),alpha10.str(8))


if __name__ == "__main__":
    verify()
