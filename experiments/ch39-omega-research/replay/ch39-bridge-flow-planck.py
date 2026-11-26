"""Validated two-loop gauge / one-loop Yukawa trajectories from the shell-matched point."""
from pathlib import Path
from runpy import run_path

from flint import acb, arb, arb_series, ctx
from sympy import fraction


def verify():
    ctx.prec = 160
    root = Path(__file__).parent
    TEN_ONE, TEN_SELF, TEN_CROSS = run_path(str(root/"ch39-two-loop-product-beta.py"))["target_inverse_coefficients"]()
    graph = run_path(str(root/"ch39-weighted-arc-zero-phase.py"))
    a, _, fit = run_path(str(root/"ch39-shell-zero-mass-hierarchy.py"))["parameters"]()
    num, den = fraction(graph["response"](3))
    def epsilon(j):
        zz = acb.zeta_zero(j)
        u = 1-1/(2*(zz.imag**2+a("1/4")))
        return graph["evaluated"](num,u)/graph["evaluated"](den,u)
    e = epsilon(11485).union(epsilon(11490))
    low,v,log_x,g5,g10 = fit(e,a(1))
    yref = (log_x.exp()*e)**2
    state = [1/g5**2, 1/g10**2, (2*log_x).exp(), arb(1)]
    pi = arb.pi()
    unit = 1/(8*pi**2)
    def rhs(state):
        g,h,x,z = state
        y = yref*z
        ca = a("54/5")/g+a("297/20")/h
        return [unit*(a("25/6")-(a("5921/15")/g+594/h-45*(x+y))*unit/2),
                unit*(TEN_ONE-(TEN_CROSS/g-TEN_SELF/h-15*(x+y))*unit/2),
                unit*x*(13*x+10*y-ca),unit*z*(10*x+13*y-ca)]

    def jet(point, degree):
        coeffs = [[value] for value in point]
        for k in range(degree):
            slope = rhs([arb_series(seq,prec=k+1) for seq in coeffs])
            for seq, derivative in zip(coeffs,slope):
                seq.append(derivative[k]/(k+1))
        return coeffs

    step = a("1/128")
    margin = (-step).union(step)
    planck_t = (a("1.2209e19")/v).log()
    crossed = {a_str: [] for a_str in ("1/4","1")}
    for i in range(1630):
        tube = [q+margin for q in state]
        deriv = rhs(tube)
        assert all(abs(f)<1 for f in deriv)
        assert deriv[0]<0 and deriv[1]>0
        if i*step <= planck_t.upper() and (i+1)*step >= planck_t.lower():
            # Re-evaluate the actual short interval, inside the same tube.
            dt = planck_t-i*step
            box = [s+dt*f for s,f in zip(state,deriv)]
            al, at = 1/(4*pi*box[0]),1/(4*pi*box[1])
            assert a(".243") < al < a(".244") and a(".0174") < at < a(".0176")
            assert a(".421") < box[2].sqrt() < a(".423")
            print("Planck:",al.str(10),at.str(10),"x",box[2].sqrt().str(10))
        coeff = jet(state,4)
        top = jet(tube,5)
        new = []
        for c,b in zip(coeff,top):
            val = b[-1]
            for ck in reversed(c):
                val = ck+step*val
            new.append(val)
        for text, times in crossed.items():
            target = 1/(4*pi*a(text))
            if state[0].upper() >= target.lower() and new[0].lower() <= target.upper():
                times.append(i)
        state = new
        if state[0] < 1/(4*pi):
            break
    else:
        raise AssertionError("strong-coupling marker not reached")
    for label, times in crossed.items():
        first,last = min(times)*step,(max(times)+1)*step
        print("alpha5 =",label,"t",first.str(10),last.str(10),
              "energy /GeV", (v*first.exp()).str(8),(v*last.exp()).str(8))
    assert (min(crossed["1/4"])*step > a("9.8") and (max(crossed["1/4"])+1)*step < a("9.9"))
    assert (min(crossed["1"])*step > a("12.0") and (max(crossed["1"])+1)*step < a("12.1"))


if __name__ == "__main__":
    verify()
