"""Sign and minimum scale of the heavy-neutrino quadratic threshold in the shell model."""
from pathlib import Path
from runpy import run_path

from flint import arb

finite = run_path(str(Path(__file__).with_name("ch39-shell-gw-finite-mass-eigenvalue.py")))
q = finite["q"]


def certify():
    z = q("-.009").union(q("-.001"))
    aux_light, _, error, _ = finite["face"](z, 5)
    before = (aux_light+error).upper()*q("1e24")
    Q = q("4.352e13").union(q("4.473e13"))
    vev = q("246.22")
    # Smallest possible largest mass across both NuFIT 6.0 IO and NO panels.
    max_floor = q(".002421").sqrt()*q("1e-9")
    # Rank-one step at the lightest sterile threshold, D=d I and d²=Q*mmax.
    first_cost = Q.lower()**3*max_floor/(4*arb.pi()**2*vev**2)
    high_light = (before-first_cost).upper()
    # For arbitrary Yukawa columns the nuclear norm bounds sum m_i
    # by sum ||D_i||²/M_i, for either ordering.
    sum_floor = (q(".0000692").sqrt()+q(".002421").sqrt())*q("1e-9")
    all_cost = Q.lower()**3*sum_floor/(4*arb.pi()**2*vev**2)
    print("pre-threshold light Rayleigh upper, minimum first step, post-threshold upper in GeV²",
          *(v.str(8) for v in (before, first_cost.lower(), high_light)))
    print("sum of all same-sign steps above", all_cost.lower().str(8))
    assert high_light < q("-1.40e24") and all_cost > q("1.9e24")


if __name__ == "__main__":
    certify()
