"""Sensitivity of the first matching threshold to the displayed gauge precision."""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx

counts = run_path(str(Path(__file__).with_name("ch39-rh-clock-solution-count.py")))
fixed = counts["fixed"]


def main():
    ctx.prec, ctx.cap = 1536, 3
    base = fixed["inputs"]()
    center = list(map(arb, (".467453", ".63811", "1.0583")))
    half = list(map(arb, (".0000005", ".000005", ".00005")))
    raw = [8*arb.pi()**2/((g-e).union(g+e))**2 for g, e in zip(center, half)]
    d, e, h = raw
    kbox = 5*d-3*e-2*h
    fbox = arb(109)/15*(e-h)-arb(23)/6*(d-e)
    # On the entire candidate critical slab, the lightest roots have
    # elasticity near one. Differentiation here holds Q = sqrt(C_1) fixed.
    unit = arb(273_000_000_000).union(arb(279_000_000_000))
    slab = counts["at_unit"]((*base[:4], kbox, fbox), unit)
    star = counts["at_first_pole"](slab)
    c1 = slab[0][0]
    mf2 = c1*(-2*star).exp()
    ratio = (125**2+6*slab[2]*mf2)/c1
    assert ratio < 1 and star > arb(24).log()/2
    # For x < C1, x = (h+L*p)/(1+L*sum(rho/(C-x))).
    floor = (125**2+min(x for x,n in slab[3])*mf2*slab[2])/(1+6*mf2/c1/(1-ratio))
    assert floor > 1000**2
    emin = (1-125**2/floor)/(1+6*mf2/c1/(1-ratio)**2)
    assert (22*23*emin-448)/16 > 3
    assert (22*23-448)/16 < 4
    # This contact is the first possible crossing for every gauge input in
    # the box. Before the decreasing branch the pure modes alone suffice.
    q0 = 1000*(kbox/44).exp()
    pure = 25*sum(fixed["positive"]((q0/c.sqrt()).log()) for c in slab[0])
    assert 22*pure-15*fbox > 0
    late = (1-arb(1)/64)/(1+6*mf2/c1/(1-ratio)**2)
    assert 22*23*((1+arb(4)/11)*late-arb(4)/11)-448 > 0
    assert slab[2]-125**2*sum(r/(c-125**2) for r,c in zip(slab[1], slab[0])) > 0
    # Thus dR/dK in (3,4), giving dR/dd_i signs (+,-,+).
    assert 5*3+15*arb(23)/6 > 0
    assert -3*3-15*(arb(109)/15+arb(23)/6) < 0
    assert -2*4+109 > 0
    for direction, lo in [(1, 273_586_101_717), (-1, 277_996_110_578)]:
        g = [v+direction*s*r for v, r, s in zip(center, half, (1,-1,1))]
        d, e, h = [8*arb.pi()**2/x**2 for x in g]
        data = (*base[:4], 5*d-3*e-2*h, arb(109)/15*(e-h)-arb(23)/6*(d-e))
        assert counts["minimum_residual"](counts["at_unit"](data, lo)) > 0
        assert counts["minimum_residual"](counts["at_unit"](data, lo+1)) < 0
        print("Gauge input corner brackets the contact unit in GeV:", lo, lo+1)
    # At the fixed clock, certify the three sign brackets throughout the full
    # display box: residual is increasing in K at fixed w and decreasing in F.
    unit_data = (*base[:4], kbox, fbox)
    fixed["monotonicity"](unit_data)
    for left, right, sign, wanted in [(".30", ".33", 1, 119), ("1.63", "1.66", -1, 71),
                                      ("8.15", "8.24", 1, 23)]:
        for w, expected_sign in [(arb(left), -sign), (arb(right), sign)]:
            k = kbox.lower() if expected_sign > 0 else kbox.upper()
            f = fbox.upper() if expected_sign > 0 else fbox.lower()
            assert expected_sign*fixed["endpoint"](w, (*base[:4], k, f))[-1] > 0
        qlo, mflo, _, _ = fixed["scales"](arb(right), kbox.lower())
        qhi, mfhi, _, _ = fixed["scales"](arb(left), kbox.upper())
        heavy = fixed["endpoint"](arb(left), (*base[:4], kbox.upper(), fbox))[-2]
        light = fixed["endpoint"](arb(right), (*base[:4], kbox.lower(), fbox))[-2]
        opened = 0
        assert mflo*arb(8).sqrt() > 1000
        for (up,n), (down,nn) in zip(heavy,light):
            assert n == nn and down > 1000
            if up < qlo:
                opened += n
            else:
                assert qhi < down
        assert opened == wanted
        print("Whole input box: w in", left, right, "active", opened,
              "Q range", qlo.str(9), qhi.str(9))


if __name__ == "__main__":
    main()
