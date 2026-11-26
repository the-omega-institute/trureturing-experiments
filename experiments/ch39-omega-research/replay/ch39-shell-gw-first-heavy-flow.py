"""Validated single-doublet matrix flow up to the first heavy threshold."""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat

folder = Path(__file__).resolve().parent
fit = run_path(str(folder/"ch39-shell-gw-orbit-strength.py"))
full = run_path(str(folder/"ch39-shell-gw-full-gauge-flavor.py"))
adj = fit["adj"]
loop = 16*arb.pi()**2
start = [arb(s)**2 for s in (".467453", ".63811", "1.0583")]
beta = [arb(41)/10, -arb(19)/6, arb(-7)]
C = [[arb(17)/20,arb(9)/4,8], [arb(1)/4,arb(9)/4,8], [arb(9)/4,arb(9)/4,0]]
eye = acb_mat([[int(i == j) for j in range(3)] for i in range(3)])
weights = [3, 3, 1]


def masses(index):
    model = full["model"]
    if index == 43:
        sing, family, _, _ = model["m"]["transported_spectra"]()
    else:
        sing, family, _ = model["scalar_models"]()
    levels = full["all_levels"](sing, family)
    inverse = [4*arb.pi()/s for s in start]
    slope = arb_mat([[beta[0]-beta[1], -arb(1)/15], [beta[1]-beta[2], arb(1)/6]])
    t, E = slope.solve(2*arb.pi()*arb_mat([[inverse[0]-inverse[1]], [inverse[1]-inverse[2]]])).entries()
    Q = 1000*t.exp()
    def offset(scale):
        return sum((Q/(125**2+scale**2*d).sqrt()).log() for d in levels)-E
    lo, hi = arb("1e12"), arb("1.2e12")
    assert offset(lo) > 0 and offset(hi) < 0
    for _ in range(18):
        mid = (lo+hi)/2
        if offset(mid)>0:
            lo = mid
        else:
            assert offset(mid)<0
            hi = mid
    result = [(125**2+lo.union(hi)**2*d).sqrt() for d in levels]
    assert 1000 < result[0] and result[-1] < Q
    return Q, result


def gauge(x, parameters):
    start, beta = parameters
    return [g/(1-2*b*g*x) for g,b in zip(start,beta)]


def jets(initial, x, order, parameters):
    """Normalized Taylor coefficients for the three Yukawa matrices in x=log(mu/TeV)/16pi²."""
    coeff = [initial]
    hats = []
    g = gauge(x, parameters)
    beta = parameters[1]
    for k in range(order):
        left = [sum((coeff[r][s]*adj(coeff[k-r][s]) for r in range(k+1)), acb_mat(3,3))
                for s in range(3)]
        trace = sum(w*sum(M[i,i].real for i in range(3)) for w,M in zip(weights,left))
        gk = [v*(2*b*v)**k for v,b in zip(g,beta)]
        shifts = [trace-sum(a*v for a,v in zip(row,gk)) for row in C]
        diff = (left[0]-left[1])*arb("1.5")
        hats.append([shifts[0]*eye+diff, shifts[1]*eye-diff, shifts[2]*eye+arb("1.5")*left[2]])
        coeff.append([sum((hats[r][s]*coeff[k-r][s] for r in range(k+1)),acb_mat(3,3))/(k+1)
                      for s in range(3)])
    return coeff


def advance(initial,x,h,parameters):
    # Picard enclosure by strict first-exit test, followed by a Taylor remainder
    # evaluated over the entire time tube. Real and imaginary components are independent.
    pad = arb(0, ".02")
    box = acb_mat([[acb(pad,pad) for j in range(3)] for i in range(3)])
    tube = [m+box for m in initial]
    span = arb(0).union(h)
    speed = jets(tube, x+span, 1,parameters)[1]
    assert all(abs(h*z)<arb(".0199") for m in speed for z in m.entries())
    tube = [a+span*b for a,b in zip(initial, speed)]
    coeff = jets(initial,x,4,parameters)
    final_coeff = jets(tube,x+span,5,parameters)[5]
    result = []
    for s in range(3):
        value = final_coeff[s]
        for term in reversed(coeff):
            value = h*value+term[s]
        result.append(value)
    return result


def inputs():
    V, u, d = fit["orbit"]["reference"]()
    lepton = acb_mat(3,3)
    for i,z in enumerate((".0000028227",".00059465",".010101")):
        lepton[i,i] = arb(z)
    return [u*arb(".8616"), V*d*arb(".01401"), lepton]


def evolve(end, initial=None, parameters=(start,beta)):
    h = arb(1)/1024
    steps = int(float((end/h).mid()))
    assert steps*h < end < (steps+1)*h
    data = inputs() if initial is None else initial
    for i in range(steps):
        data = advance(data,i*h,h,parameters)
    return advance(data,steps*h,end-steps*h,parameters)


def at_threshold(index):
    Q, heavy = masses(index)
    scale = (heavy[0]/1000).log()/loop
    up, down, lepton = evolve(scale)
    F, D = fit["readouts"](index)
    lo, hi, field = fit["limits"](index,F,D,up,down)
    a,b = (("5.17464","5.17551") if index == 43 else ("5.09736","5.09836"))
    assert arb(a) < lo and hi < arb(b)
    # Compare the very same frozen color-singlet principal block for this witness.
    assert fit["count"](fit["rows"](F*adj(F)),arb(".52")**2) == 9
    assert 3*(arb(".52")*hi)**2/(16*arb.pi()) < arb(".437")
    print("zero",index,"Q",Q.str(8),"first heavy GeV",heavy[0].str(12),
          "norms",*[fit["hs"](m,m).real.sqrt().str(12) for m in (up, down)],
          "common strength lower and witness",lo.str(12),hi.str(12))
    return heavy, up, down, lepton, F, D, lo, hi, field


if __name__ == "__main__":
    for n in (43,195):
        at_threshold(n)
