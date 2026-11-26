"""Test the leading quartic feasibility at large arc inertia in the geometric zero window."""
from fractions import Fraction
from pathlib import Path
import runpy

from flint import arb

root = Path(__file__).resolve().parent
model = runpy.run_path(str(root / "ch39-quaternion-scalar-levels.py"))
pencil = model["m"]["m"]
geometry = pencil["real_energy_matrix"].__globals__["src"]
geometry["rho"] = arb.pi()**2/geometry["theta"](804)
original = pencil["real_energy_matrix"]
edge_types = [d for _, _, d, _ in geometry["edges"]]
angles = {d: geometry["theta"](d)/geometry["theta"](804) for d in set(edge_types)}
q = lambda n: model["q"](Fraction(n))


def eigen_count(energy, vertex_mass, direction):
    # E/d is the spectral variable. The sine and cosine terms have no d left.
    def matrix(s):
        H = original(s)
        clock = [arb.pi()*(s*angles[d]).sqrt() for d in (724, 804)]
        shift = q(".04")*max(abs(k/k.sin()).upper() for k in clock)
        for i in range(H.nrows()):
            H[i, i] += (1-vertex_mass)*s + direction*shift
        return H

    pencil["real_energy_matrix"] = matrix
    dn = sum(int((energy*angles[d]).sqrt().floor().unique_fmpz()) for d in edge_types)
    return 2*dn + model["count"](energy)


if __name__ == "__main__":
    # Uniform lower mass at d <= 1/2 via monotonicity of E/(d*P+C).
    dmax = q("0.5")
    h = q("1e-6")
    trial = q(".0585")
    assert eigen_count(trial, dmax, -1) == 0
    coeff = arb.pi()**2*h**2/6
    first_lower = trial/(1+coeff*(trial+10/dmax))
    # Upper ratio via 128 spectral counts for the continuum at zero vertex mass.
    ratio = q("1.44")
    edges = [q(".20001")+(ratio-q(".20001"))*j/128 for j in range(129)]
    budget = sum(eigen_count(a*a, q(0), 1)*(b/a).log() for a,b in zip(edges,edges[1:]))
    assert budget > q("49.2")
    # Interpolation with controlled endpoint traces, even as d -> 0.
    assert min(arb.pi()**2*z for z in angles.values()) > 1
    trace_bound = 80 + q(".1")*ratio**2
    retained = (1-h**2*ratio**2)**2 - 9*arb.pi()**2*h*trace_bound
    upper_ratio = ratio/retained.sqrt()
    width = (upper_ratio/first_lower.sqrt()).log()
    lh = q(".02335")
    lam_upper = q("-.0298109")
    target = 2*arb.pi()**2*(lh-lam_upper)/3
    eta = 1+4*lh/(8*arb.pi()-6*lh)
    U = 16*arb.pi()*lh
    required = 6*target*eta/U
    assert width < required
    # In case the stability/contact test is made after adding a common b,
    # its positive part must cover this Hilbert-Schmidt norm deficit.
    # Each mass level has at most six initially negative complex directions.
    trace_shift = ((6*target/width).sqrt() - (U/eta).sqrt())/q(6*308).sqrt()
    assert trace_shift > 0
    for name, val in (("continuum exposure lower",budget), ("finite first lower",first_lower),
                      ("finite ratio upper",upper_ratio), ("finite width upper",width),
                      ("required width lower",required), ("necessary common b",trace_shift)):
        print(name,val.lower().str(13),val.upper().str(13))
