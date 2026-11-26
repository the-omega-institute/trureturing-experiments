"""Five shell steps are optimal within the positive sum of filtered mixed densities."""
from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat
import numpy as np

finite = run_path(str(Path(__file__).with_name("ch39-shell-gw-finite-mass-eigenvalue.py")))
q, base = finite["q"], finite["base"]


def family_levels():
    full = base["flow"]["full"]
    single, family, _, _ = full["model"]["m"]["transported_spectra"]()
    rows, count = full["rows"], full["count"]
    ref = np.linalg.eigvalsh(full["model"]["array"](single))[0]
    delta = q("1e-9")
    lower, upper = arb(float(ref))-delta, arb(float(ref))+delta
    assert count(rows(single), lower) == 0 and count(rows(single), upper) == 1
    ground = lower.union(upper)
    levels = []
    for k, val in enumerate(np.linalg.eigvalsh(full["model"]["array"](family))):
        a, b = arb(float(val))-delta, arb(float(val))+delta
        assert count(rows(family), a) == k and count(rows(family), b) == k+1
        levels.append(a.union(b)-ground+arb(0, "12e-8"))
    return levels


def parts(x, d, e, f):
    product = (x-d)*(x-e)
    return 1-product/((f-d)*(f-e)), product*(x-f), (x-4)/3


def basis(x, third):
    y = (x-4)/3
    return [(x-third)*p for p in (1, y, 2*y*y-1, 4*y**3-3*y)]


def outer(a, b):
    return arb_mat([[x*y for y in b] for x in a])


def positive(matrix):
    trial = np.array([[float(matrix[i, j].mid()) for j in range(4)] for i in range(4)])
    _, V = np.linalg.eigh(trial)
    V = arb_mat([[arb(float(v)) for v in row] for row in V])
    matrix = V.transpose()*matrix*V
    assert not V.det().contains(0)
    for n in range(1, 5):
        minor = arb_mat([[matrix[i, j] for j in range(n)] for i in range(n)])
        assert minor.det() > 0


def certify():
    levels = family_levels()
    d, e, f = levels[:3]
    assert d < e < f
    gram = arb_mat(4, 4)
    weights, fifth = [], []
    for x in levels:
        p, a, y = parts(x, d, e, f)
        weight = (p+a*(q(".25946613")-q(".12032577")*y)).mid()**2
        weights.append(weight)
        vec = basis(x, f)
        gram += weight*outer(vec, vec)
        fifth.append(p+a*(q(".33941785")-q(".24571575")*y+q(".05293755")*(2*y**2-1)))
    # Divided difference of t(x)t(x)^T on [d,e] by the mean value theorem.
    # The exact divided difference enters the quadratic-form constraint.
    x = d.union(e)
    vec = basis(x, f)
    y = (x-4)/3
    deriv = [v+(x-f)*dv for v, dv in zip((1,y,2*y*y-1,4*y**3-3*y),
                                          (0,q(1,3),4*y/3,4*y*y-1))]
    diff = outer(deriv, vec)+outer(vec, deriv)
    positive(gram-100*outer(basis(d, f), basis(d, f))-50*diff)
    normw = sum(v*v for v in weights).sqrt()
    lower = q(6).sqrt()*(2*q(".61703")*100-q(".01150001")*sum(weights))/normw
    assert lower > q("28.1") and lower > 8*arb.pi()

    events = []
    z = q("-.00675").union(q("-.0035"))
    base["check"](q(81, 977)**5, mass_path=events, bulk_portal=z)
    quart = next(row[3] for row in events if row[:2] == ("threshold", 3))
    ao, bo, bn = quart["a", (0, 1)], quart["b", (0, 1)], quart["b", (0, 2)]
    assert ao-z > q(".61703") and 2*z+bn > q("-.01150001")
    # The first two gains are exactly one and the third is exactly zero.
    sup = max(abs(v).upper() for v in fifth[3:])
    assert sup < q("1")
    fourth = 2+sum(v**4 for v in fifth[3:])
    # ||F^*P_v F|| on the family-pair output is obtained from its Gram
    # matrix in v. Row sums equal Tr(P_v F^4) by site transitivity.
    scale = (fourth/6).sqrt()
    cross = 6*(q(8).sqrt()*(2*abs(z)+abs(bn))+(2*abs(ao-z)+abs(bo-bn))*scale)
    single, family = 36*quart["l", 0], 204*quart["l", 2]
    upper = (single+family+((single-family)**2+4*cross**2).sqrt())/2
    mixed = 6*(abs(z)+3*abs(bn)+abs(ao-z)+3*abs(bo-bn))
    upper = arb.max(upper, mixed)
    assert upper < q("9.6")
    print("up to four steps: contact norm exceeds", lower.lower().str(9))
    print("five-step certificate: higher-family gain below", sup.str(8),
          "fourth spectral moment below", fourth.upper().str(8),
          "contact norm below", upper.upper().str(8))


if __name__ == "__main__":
    certify()
