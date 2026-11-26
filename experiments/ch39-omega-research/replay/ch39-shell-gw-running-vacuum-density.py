"""Vacuum-density matching and full one-loop scalar mass transport across the shell."""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat

slope = run_path(str(Path(__file__).with_name("ch39-shell-gw-finite-portal-slope.py")))
q, base, signed, unit = (slope[k] for k in ("q","base","signed","unit"))


def certify():
    events = []
    z = q("-.00675").union(q("-.0035"))
    T, B, A, margin, raw, time = base["check"](q(81, 977)**5, mass_path=events, bulk_portal=z)
    Y = (12*T).sqrt()
    heavy = q(12).sqrt()*raw.sqrt()*q(81, 977)**5*(6*T*time).exp()
    zrate = 12*A+3*q(".020")+4*Y**2
    zsource = (8*Y*B+96*Y**3)*heavy
    assert slope["propagate"](q(0),zrate,zsource,time) < min(6*margin, q(".02"))
    low = q("-1e6").union(q("1e6"))/unit
    # The helper only uses the individual and sum mass vectors in its first pair;
    # the gradient pair receives the same boundary derivative as in the portal study.
    state = ([low], [low], [q(0)], [q(0)])
    direction = {}
    elapsed = mass_error = jumps = run_lo = run_hi = q(0)
    for op, *args in events:
        if op == "threshold":
            slope["add_at_threshold"](*args, direction, state)
            scale = args[-1]
            jumps += 3*(scale**2/unit)**2/(32*arb.pi()**2)
            continue
        groups, quartic, gs, ym, dt = args
        gen = slope["coefficients"](groups, quartic, gs, ym)
        vector = arb_mat([[v] for part in state[:2] for v in part])
        n = len(groups)
        tube = (signed["advance_matrix"](gen, q(0).union(dt))*vector).entries()[n:]
        end = (signed["advance_matrix"](gen, dt)*vector).entries()
        state[0][:], state[1][:] = end[:n], end[n:]
        # Square the lower and upper absolute endpoints separately; includes zero.
        lower = 2*sum(num*tube[i].abs_lower()**2 for i,num in groups.items()).sqrt()
        upper = 2*sum(num*tube[i].abs_upper()**2 for i,num in groups.items()).sqrt()
        if n == 3:
            elapsed += dt.upper()
            quartic_error = slope["propagate"](q(0), zrate, zsource, arb.min(elapsed, time.upper()))
            mass_error = slope["propagate"](mass_error, 2*A+quartic_error+2*Y**2,
                                            (quartic_error+4*Y*heavy)*upper, dt)
        norm_lo = arb.max(q(0), (lower-mass_error).lower())
        norm_hi = (upper+mass_error).upper()
        # Normalized RG time; origin of unbroken phase, fermions and SM vectors massless.
        run_lo += arb.max(q(0), dt.lower())*norm_lo**2/2
        run_hi += dt.upper()*norm_hi**2/2
    running = run_lo.lower().union(run_hi.upper())
    total = jumps+running
    p = q("2.435e18")
    target = 8*arb.pi()*q("1e-122")*p**4
    crossing = q("4.352e13").union(q("4.473e13"))
    extra = 3*17*crossing**4/(64*arb.pi()**2)
    print("finite thresholds and integrated scalar running in 1e48 GeV^4",
          [(float(v.lower()), float(v.upper())) for v in (jumps, running, total)])
    print("17 charged thresholds in GeV^4", float(extra.lower()), float(extra.upper()))
    print("fractional vacuum-density width for ±target, shell only and with extra fields below",
          (2*target/(total.lower()*unit**2)).upper().str(8),
          (2*target/(total.lower()*unit**2+extra.lower())).upper().str(8))
    assert running > 0 and total > 50 and extra > q("2.8e53")


if __name__ == "__main__":
    certify()
