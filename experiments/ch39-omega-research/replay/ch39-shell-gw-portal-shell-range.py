"""Third-threshold local quartic obstruction and two-step shell-filter realization."""
from pathlib import Path
from runpy import run_path

finite = run_path(str(Path(__file__).with_name("ch39-shell-gw-finite-mass-eigenvalue.py")))
q = finite["q"]
base = finite["base"]


def check():
    full = base["flow"]["full"]
    single, family, _, _ = full["model"]["m"]["transported_spectra"]()
    rows, count = full["rows"], full["count"]
    shift = q("12e-8")
    bands = [(".0158803", ".0158810"), (".0158828", ".0158835"), ("1.08320", "1.08322")]
    for number, (left, right) in enumerate(bands):
        lo, hi = q(left), q(right)
        assert count(rows(family), lo+shift) == number
        assert count(rows(family), hi-shift) == number+1
        assert count(rows(single), lo) == count(rows(single), hi) == 1
    denominator = (q(bands[2][0])-q(bands[0][1]))*(q(bands[2][0])-q(bands[1][1]))
    assert denominator > q("1.139")
    events = []
    z = q("-.00675").union(q("-.0035"))
    base["check"](q(81, 977)**5, mass_path=events, bulk_portal=z)
    third = next(row for row in events if row[:2] == ("threshold", 3))
    quartic = third[3]
    a_low, b_low = (quartic[t, (0, 1)] for t in ("a", "b"))
    b_high = quartic["b", (0, 2)]
    l0, l3 = quartic["l", 0], quartic["l", 2]
    assert a_low > q(".61353") and b_low > b_high
    assert (z+b_high).upper() < q("-.00149")
    assert 6*(a_low-z) > q("3.7021")
    # Copositivity of the two unfiltered local densities with
    # coefficients 6*l0, 6*l3 and 6*z.
    margin = z+2*(l0*l3).sqrt()
    assert margin > q(".0362")
    print("a,b for the old species", a_low.str(6), b_low.str(6),
          "normalized local directional diameter above", (6*(a_low-z)).lower().str(8))
    print("positive local-potential margin", margin.lower().str(8),
          "denominator for the degree-two filter above", denominator.lower().str(8))


if __name__ == "__main__":
    check()
