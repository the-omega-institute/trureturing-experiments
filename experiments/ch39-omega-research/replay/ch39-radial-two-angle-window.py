"""Continue the least branch over both input-angle windows and certify its extrema."""
from itertools import product
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat

fit = run_path(str(Path(__file__).with_name("ch39-radial-gain-two-angle-fit.py")))
model = fit["model"]
angles = (model["s"], model["t"])
width = model["shell"]["sine_radius"]
windows = [a+arb(0,width.upper()) for a in angles]
gain = arb(18)/25


def observables_and_slopes(w, derivatives, inputs):
    ud,ub,cd,cb,us = [w[i,j] for i,j in ((0,0),(0,2),(1,0),(1,2),(0,1))]
    u, v, k = abs(ub)**2, abs(ud)**2, abs(cb)**2
    signed = -ud*ub.conjugate()*cd.conjugate()*cb
    z = signed-v*k*u/(1-u)
    conditions = [abs(x)**2-y*y*(1-u) for x,y in zip((us,cb), inputs)]
    cols, out = [], []
    def square_gradient(x, dx):
        return 2*(x.conjugate()*dx).real
    for dw in derivatives:
        dud,dub,dcd,dcb,dus = [dw[i,j] for i,j in ((0,0),(0,2),(1,0),(1,2),(0,1))]
        du,dv,dk = (square_gradient(x,dx) for x,dx in ((ub,dub),(ud,dud),(cb,dcb)))
        dz = signed*(dud/ud+(dub/ub).conjugate()+(dcd/cd).conjugate()+dcb/cb)
        dz -= (dv*k+v*dk)*u/(1-u)+v*k*du/(1-u)**2
        cols.append([square_gradient(us,dus)+inputs[0]**2*du,dk+inputs[1]**2*du])
        out.append([du/(2*abs(ub)), (dz/z).imag])
    return conditions, arb_mat(cols).transpose(), arb_mat(out).transpose()


def regular_box(strength, phase, response):
    # Exclude boxes whose model cannot produce both input angles.
    w, derivatives = model["mixing_matrix"]((gain,strength,phase),True,response)
    conditions, jac, outputs = observables_and_slopes(w, derivatives[1:], windows)
    if any(0 not in v for v in conditions):
        return True
    determinant = jac[0,0]*jac[1,1] - jac[0,1]*jac[1,0]
    if not determinant.is_finite() or 0 in determinant:
        return False
    adjugate = arb_mat([[jac[1,1], -jac[0,1]], [-jac[1,0], jac[0,0]]])
    slope_signs = outputs*adjugate/determinant
    # The two suppressed factors 2*s_i*(1-|Vub|^2) are strictly positive.
    return (slope_signs[0,0]>0 and slope_signs[0,1]>0
            and slope_signs[1,0]>0 and slope_signs[1,1]<0)


def cover(response):
    # First reduce the strength interval using both phase quadratures.
    pending = [(arb(0),arb("2.6"))]
    strength_hull = phase_hull = None
    while pending:
        lo,hi = pending.pop()
        h = lo.union(hi)
        possible, phase = fit["possible_phase"](h,response,gain,windows)
        if not possible:
            continue
        if (phase is not None and h>arb("2.54") and h<arb("2.58")
                and phase>arb("3.04") and phase<arb("3.13")):
            strength_hull = h if strength_hull is None else strength_hull.union(h)
            phase_hull = phase if phase_hull is None else phase_hull.union(phase)
            continue
        assert hi-lo > arb("1e-11")
        mid = ((lo+hi)/2).mid()
        pending.extend(((lo,mid),(mid,hi)))
    count = 0
    assert strength_hull is not None and phase_hull is not None
    boxes = [(strength_hull,phase_hull)]
    while boxes:
        h, phase = boxes.pop()
        count += 1
        if regular_box(h, phase, response):
            continue
        where = 0 if h.rad() >= phase.rad() else 1
        x = (h,phase)[where]
        assert x.rad() > arb("1e-10")
        halves = (x.lower().union(x.mid()), x.mid().union(x.upper()))
        boxes.extend((part,phase) if where==0 else (h,part) for part in halves)
    return count


def corner(inputs, response):
    point = [arb("2.56051"),arb("3.08330")]
    def system(box):
        w,d = model["mixing_matrix"]((gain,*box),True,response)
        f,j,_ = observables_and_slopes(w,d[1:],inputs)
        return arb_mat([[v] for v in f]), j
    for _ in range(4):
        f,j = system(point)
        corrections = j.solve(f)
        point = [(v-corrections[i,0]).mid() for i,v in enumerate(point)]
    r = arb("1e-8")
    box = [v+arb(0,r) for v in point]
    f,j = system(point)
    inverse = j.inv()
    rem = fit["eye"]-inverse*system(box)[1]
    lip = max(sum(abs(rem[i,j]).upper() for j in range(2)) for i in range(2))
    assert max(abs(v).upper() for v in (inverse*f).entries())+lip*r < r
    w = model["mixing_matrix"]((gain,*box),response=response)
    ud,ub,cd,cb = [w[i,j] for i,j in ((0,0),(0,2),(1,0),(1,2))]
    z = -ud*ub.conjugate()*cd.conjugate()*cb-abs(ud*cb*ub)**2/(1-abs(ub)**2)
    return abs(ub), z.arg()


def verify(response):
    outer, _ = fit["isolate"](("2.56051","3.0833"),response,gain)
    fit["prove_least"](response,gain,outer)
    count = cover(response)
    choices = dict()
    for signs in product((-1,1),repeat=2):
        choices[signs] = corner([a+sign*width for a,sign in zip(angles,signs)],response)
    minimum_u, maximum_u = choices[-1,-1][0], choices[1,1][0]
    minimum_phase, maximum_phase = choices[-1,1][1], choices[1,-1][1]
    assert arb(".003746")<minimum_u and maximum_u<arb(".003780")
    assert arb("1.1436")<minimum_phase and maximum_phase<arb("1.1452")
    print("certified local boxes",count,"s13 extrema",minimum_u.str(8),maximum_u.str(8))
    print("least branch over both full input windows; CP extrema",minimum_phase.str(8),maximum_phase.str(8))


if __name__ == "__main__":
    verify(fit["circle"]["two_ticks"]["verify"]())
