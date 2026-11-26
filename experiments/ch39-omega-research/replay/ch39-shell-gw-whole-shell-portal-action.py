"""Contact unitarity on the full shell for a polynomial light–family portal."""
from pathlib import Path
from runpy import run_path

from flint import arb

finite = run_path(str(Path(__file__).with_name("ch39-shell-gw-finite-mass-eigenvalue.py")))
q, base = finite["q"], finite["base"]


def chebyshev(n, x):
    previous, current = q(1), x
    for _ in range(1, n):
        previous, current = current, 2*x*current-previous
    return current


def checks():
    full = base["flow"]["full"]
    single, family, _, _ = full["model"]["m"]["transported_spectra"]()
    rows, count = full["rows"], full["count"]
    # Ordered eigenvalues of the shifted family block on the entire angle square.
    shift = q("12e-8")
    lo0, hi0 = q(".01375"), q(".01377")
    assert count(rows(single), lo0) == 0 and count(rows(single), hi0) == 1
    for k, (lo, hi) in enumerate(((".0158803", ".0158810"), (".0158828", ".0158835"),
                                 ("1.08320", "1.08322"))):
        assert count(rows(family), q(lo)+shift) == k
        assert count(rows(family), q(hi)-shift) == k+1
    assert count(rows(family), q(6)+hi0+shift) < 48
    assert count(rows(family), q(7)+lo0-shift) == 48
    low = q(".0158803").union(q(".0158835")) - lo0.union(hi0)
    third = q("1.08320").union(q("1.08322")) - lo0.union(hi0)
    max_lower = q(6)
    # Convexity of the norm of an affine vector below a root shared by
    # all the filters bounds their gain beyond that root.
    minimum_gain = (max_lower-third)*(third+max_lower-2*low)/(third-low)**2
    too_large = q("-.00675")+(q(".61353")+q(".0035"))*minimum_gain**2
    assert too_large > 560

    events = []
    z = q("-.00675").union(q("-.0035"))
    base["check"](q(81, 977)**5, mass_path=events, bulk_portal=z)
    quartic = next(row[3] for row in events if row[:2] == ("threshold", 3))
    l0, l3 = quartic["l", 0], quartic["l", 2]
    ao, bo = quartic["a", (0, 1)], quartic["b", (0, 1)]
    bn = quartic["b", (0, 2)]
    # Interpolate 1 / ((d3-x)*T_n((4-x)/3)) at the two low nodes.
    # The divided difference is bounded by the derivative on their hull.
    n = 9
    y = (4-low)/3
    # Hyperbolic evaluation avoids dependency in the low-frequency tube.
    angle = y.acosh()
    T = (n*angle).cosh()
    g = 1/((third-low)*T)
    gprime = g*(1/(third-low)+n*(n*angle).tanh()/(3*(y**2-1).sqrt()))
    assert gprime > 0
    eps = (7-third)*(g+(7-low)*gprime)
    assert eps < q(".31")
    local_gain = (q(1,3)+eps**2).sqrt()
    # Same-species symmetric two-body blocks; two light and two family
    # legs occur in every mixed term.
    light_diag = 36*l0
    family_diag = 204*l3
    base_cross = q(8).sqrt()*6*(2*abs(z)+abs(bn))
    filtered_cross = q(8).sqrt()*6*(2*abs(ao-z)+abs(bo-bn))*local_gain
    total_cross = base_cross+filtered_cross
    diagonal_bound = (light_diag+family_diag+((light_diag-family_diag)**2
                      +4*total_cross**2).sqrt())/2
    mixed_bound = 6*(abs(z)+3*abs(bn)+abs(ao-z)+3*abs(bo-bn))
    upper = arb.max(diagonal_bound, mixed_bound)
    print("any two-step common-floor portal above", too_large.lower().str(9),
          "eleven-step stop-band gain below", eps.upper().str(8))
    print("full scalar contact norm below", upper.upper().str(8),
          "normalized scalar J=0 norm below", (upper/(16*arb.pi())).upper().str(8))
    assert upper < q("18.0")
    return upper


if __name__ == "__main__":
    checks()
