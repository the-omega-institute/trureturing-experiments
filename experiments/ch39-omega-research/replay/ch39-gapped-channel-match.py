"""Universal matching bounds over positive channel spectra with a common gap floor."""
from pathlib import Path
from runpy import run_path

from flint import arb

contact = run_path(str(Path(__file__).with_name("ch39-624-clock-contact.py")))
clock, mixture, flavor = (contact[k] for k in ("clock", "mixture", "flavor"))


def verify():
    rows = contact["clock_window"]()
    _,p0,w0,reference = next(row for row in rows if row[0] == 11431)
    others = [(p0-p,w/w0,[a-b for a,b in zip(obs,reference)]) for n,p,w,obs in rows if n!=11431]
    e0 = clock["free_ground"](p0)
    def boundary(t):
        energy = clock["free_ground"](p0-t)
        a,da = clock["endpoint"](energy,1)
        b,db = clock["endpoint"](energy,3)
        return a-b, db+p0*(da-db), da-db
    _,slope,_ = boundary(arb(0))
    residue = w0/(4*slope)
    def secular(t):
        span,_,_ = boundary(t)
        return 4*span*(625-arb(1)/2)*t/w0-1+sum(w*t/(d-t) for d,w,_ in others)
    lo,hi = arb(0), (min(d.lower() for d,_,_ in others)/2).mid()
    assert secular(hi)>0
    for _ in range(26):
        mid = ((lo+hi)/2).mid()
        if secular(mid)<0:
            lo=mid
        else:
            assert secular(mid)>0
            hi=mid
    limit = clock["free_ground"](p0-hi)
    assert limit < min(clock["free_ground"](p0-d).lower() for d,_,_ in others)
    assert limit-e0 < arb(".000001")

    # For any excited channel, E-level is below this endpoint. The squared
    # resolvent norm grows monotonically below the first pole of each block.
    floor = arb(".1")
    a,da = clock["endpoint"](limit-floor,1)
    b,db = clock["endpoint"](limit-floor,3)
    high_bounds = []
    for n,p,w,obs in rows:
        m = b+p*(a-b)
        derivative = (db+p*(da-db))/(1-2*m)**2
        high_bounds.append((624*w*derivative,obs))
    high_mass = sum(w for w,_ in high_bounds)
    high_offsets = [sum(w*arb(0).union(obs[k]-c) for w,obs in high_bounds)
                    for k,c in enumerate(reference)]
    radii = [flavor["shell"]["sine_radius"]]*3+[flavor["shell"]["phase_radius"]]
    intervals = [(arb(0),hi)]
    worst_gauge = arb(0)
    least_probability = arb(1)
    while intervals:
        lo,stop = intervals.pop()
        t = lo.union(stop)
        _,ds,dd = boundary(t)
        norm = arb(1)
        low_offsets = [arb(0) for _ in reference]
        for d,w,delta in others:
            ratio = w*(1-d*dd/ds)*(t/(d-t))**2
            norm += ratio
            for k,f in enumerate(delta):
                low_offsets[k] += ratio*f
        scale = ((clock["free_ground"](p0-stop)-e0)**2/residue).upper()
        # Normalize the high-sector bounds by the full low-mode norm.
        mass_bound = scale*high_mass
        denom = 1+arb(0).union(mass_bound)
        obs = [c+(x/norm+arb(0).union(scale)*y)/denom
               for c,x,y in zip(reference,low_offsets,high_offsets)]
        vals = mixture["readouts"](obs)
        drift = obs[-1]/reference[-1]-1
        error = arb("1.005e-7")+arb(".13")*abs(drift.log1p())/2
        if error < arb(2)**-21 and all(abs(v-c)<r for v,c,r in zip(vals,flavor["shell"]["centers"],radii)):
            worst_gauge = max(worst_gauge,error.upper())
            least_probability = min(least_probability,(1/(norm*denom)).lower())
            continue
        assert stop-lo > arb("1e-12")
        mid = ((lo+stop)/2).mid()
        intervals.extend(((lo,mid),(mid,stop)))
    assert least_probability > arb(".999857") and worst_gauge < arb("4.22e-7")
    print("every 625-state uniform-ground channel, gap >=",floor.str(3),
          "11431 probability >",least_probability.str(8),"three residuals <",worst_gauge.str(8))


if __name__=="__main__":
    verify()
