"""Whole-shell onsite positive quartic extensions on the first eight active layers."""
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import acb, acb_mat, arb

finite = run_path(str(Path(__file__).with_name("ch39-shell-gw-finite-mass-eigenvalue.py")))
q, base = finite["q"], finite["base"]
factors = [q(1), q("1.02"), q("2.2")]


def lower_gram():
    full = base["flow"]["full"]
    single, family, _, _ = full["model"]["m"]["transported_spectra"]()
    rows, count = full["rows"], full["count"]
    data = rows(family)
    assert count(rows(single), q("2.29")) == 1
    assert count(data, q(".04")) == count(data, q("1.06")) == 2
    assert count(data, q("1.13")) == count(data, q("2.28")) == 8
    vals, vecs = np.linalg.eigh(full["model"]["array"](family))
    vals = [arb(float(val)) for val in vals[:8]]
    assert max(vals[:2]) < q(".016") and min(vals[2:]) > q("1.08")
    assert max(vals[2:]) < q("1.12")
    W = acb_mat([[acb(float(z.real), float(z.imag)) for z in row] for row in vecs[:, :8]])
    adj = lambda a: a.conjugate().transpose()
    eye = acb_mat([[int(i == j) for j in range(8)] for i in range(8)])
    defect = adj(W)*W-eye
    assert sum(abs(z).abs_upper()**2 for z in defect.entries()).sqrt() < q("1e-12")
    diagonal = acb_mat([[vals[i] if i == j else 0 for j in range(8)] for i in range(8)])
    residual = family*W-W*diagonal
    assert sum(abs(z).abs_upper()**2 for z in residual.entries()).sqrt() < q("1e-12")
    # Cluster projections for each group, followed by separate polar factors.
    delta = q("12e-8")
    separate = (8*q(1)).sqrt()*(delta*(1+q("1e-12"))+q("1e-12"))
    distorted = q("1e-12")+2*separate*(1+q("1e-12"))+separate**2
    # Frobenius distance to the two column blocks' polar corrections.
    error = separate+q(8).sqrt()*distorted/(1+(1-distorted).sqrt())
    gram_error = 6*error*(2*(1+q("1e-12"))+error)
    assert gram_error < q(".00004")
    for v in range(6):
        block = acb_mat(rows(W)[8*v:8*v+8])
        G = 6*adj(block)*block
        for j in range(8):
            G[j, j] -= 1/factors[1 if j < 2 else 2]
        assert count(rows(G), q(".014")) == 0


def weighted_contact(groups, quartic):
    largest = q(0)
    for i, count in groups.items():
        row = 6*factors[i]**2*abs(quartic["l", i])
        for j, size in groups.items():
            copies = size-int(i == j)
            if copies:
                ij = base["pair"](i, j)
                a = abs(quartic["a", ij])*factors[i]*factors[j]
                b = abs(quartic["b", ij])*factors[i]*factors[j]
                row += copies*(2*a+b)
                largest = arb.max(largest, (a+3*b).upper())
        largest = arb.max(largest, row.upper())
    return largest


def certify():
    lower_gram()
    events = []
    z = q("-.00675").union(q("-.0035"))
    T, B, A, margin, raw, time = base["check"](q(81, 977)**5, mass_path=events, bulk_portal=z)
    Y = (12*T).sqrt()
    heavy = q(12).sqrt()*raw.sqrt()*q(81, 977)**5*(6*T*time).exp()
    rate = 12*A+3*q(".020")+4*Y**2
    error = (8*Y*B+96*Y**3)*heavy*((rate*time).exp()-1)/rate
    assert error < min(6*margin, q(".02"))
    norm = q(0)
    for op, *args in events:
        if op == "threshold":
            k, groups, quartic, _ = args
            if k > 8:
                break
        else:
            groups, quartic, _, _, _ = args
        norm = arb.max(norm, weighted_contact(groups, quartic))
    # A separate positive potential on the orthogonal complement of each
    # evaluation image adds a disjoint two-particle block.
    actual = 6*(norm+factors[2]**2*error/2)
    assert actual < 8*arb.pi()
    print("whole-shell onsite contact norm on all eight layers below", actual.upper().str(9),
          "J=0 norm below", (actual/(16*arb.pi())).upper().str(8))


if __name__ == "__main__":
    certify()
