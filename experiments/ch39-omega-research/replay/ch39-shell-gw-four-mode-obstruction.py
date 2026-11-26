"""Certified exclusion of the four lowest charged-shell Higgs modes for quark mixing."""
from itertools import permutations
from pathlib import Path
from runpy import run_path

import numpy as np
from flint import acb, acb_mat, arb

band = run_path(str(Path(__file__).with_name("ch39-shell-gw-higgs-full-band.py")))
hs, adj, frame = (band[name] for name in ("hs", "adj", "frame"))
yuk = band["yuk"]


def orthonormalize(vectors):
    result = []
    for vector in vectors:
        v = vector
        for b in result:
            v = v - hs(b, v) * b
        n = hs(v, v).real
        assert n > 0
        result.append(v / n.sqrt())
    return result


def float_matrix(a):
    return np.array([[complex(a[i, j]) for j in range(a.ncols())]
                     for i in range(a.nrows())])


def unitary_trial(a):
    u, _, vh = np.linalg.svd(a)
    return u @ vh


def exact_frame(trials):
    cols = [acb_mat([[acb(float(z.real), float(z.imag))] for z in col])
            for col in trials.T]
    axes = orthonormalize(cols)
    return acb_mat([[col[i, 0] for col in axes] for i in range(3)])


def four_mode_basis():
    sing, family, phases, generators = (band[k] for k in ("sing", "family", "phases", "g"))
    count, rows = band["count"], band["rows"]
    fifth_lo, fifth_hi = arb("1.0868134714"), arb("1.0868134716")
    assert all(count(rows(sing), c) == 1 for c in (arb(".02"), fifth_hi))
    assert [count(rows(family), c) for c in (arb(".02"), arb("1.085"), fifth_lo, fifth_hi)] == [2, 3, 3, 4]
    singlet = [yuk["identity"](3) / arb(3).sqrt()]
    sources = [(frame(sing, arb(".02"))[0], singlet)]
    sources += [(v, generators) for v in frame(family, arb(".02"))]
    sources.append((frame(family, arb("1.085"))[2], generators))
    lzero, rzero = (yuk["zero_basis"](s) for s in (1, -1))

    def vertex(a, i):
        return acb_mat([[a[3*i+r, c] for c in range(3)] for r in range(3)])

    matrices = []
    for v, basis in sources:
        y = acb_mat(3, 3)
        for i in range(6):
            h = sum((v[i*len(basis)+j, 0]*g for j, g in enumerate(basis)), acb_mat(3, 3))
            y += phases[i] * vertex(lzero, i).transpose() * h.transpose() * vertex(rzero, i)
        matrices.append(y)
    return orthonormalize(matrices), matrices


def certify_entries():
    modes, raw = four_mode_basis()
    y0, y1 = map(float_matrix, raw[:2])
    _, eigenvectors = np.linalg.eig(y1 @ np.linalg.inv(y0))
    left = unitary_trial(eigenvectors)
    right = unitary_trial(y0.conj().T @ left)
    L, R = exact_frame(left), exact_frame(right)
    modes = [adj(L) * y * R for y in modes]
    for i, j, _ in permutations(range(3)):
        assert sum(abs(b[i, j])**2 for b in modes[:3]).sqrt() < arb(".000415")
    # The eigensolver may order the three axes in either cyclic orientation.
    forward = [(0, 1), (1, 2), (2, 0)]
    if abs(modes[3][0, 1]) < arb(".2"):
        forward = [(j, i) for i, j in forward]
    for i, j in forward:
        assert arb(".563") < abs(modes[3][i, j]) < arb(".565")
        assert abs(modes[3][j, i]) < arb(".125")


def heavy_line_radius(r, delta):
    """Bound the projector distance of a leading left singular line from one fixed axis."""
    a, A, B = arb(".563"), arb(".565"), arb(".125")
    # The two cyclic products of entries of every rank-one matrix agree.
    t = 2 * (r + delta) / (a - B)
    for _ in range(3):
        off = (A+B)*t + 2*delta
        assert off < 1
        # Minor through the largest diagonal entry and the disjoint forward edge.
        t = ((r+(B*t+delta)**2)/(1-off)+delta)/a
    off = (A+B)*t + 2*delta
    coupling = ((A*t+delta)**2+(B*t+delta)**2).sqrt()
    assert r + 2*off < 1
    return coupling/(1-r-2*off)


if __name__ == "__main__":
    certify_entries()
    tol = arb("1e-6")
    ratio_errors = (1+tol)/(1-tol)
    ru, rd = arb(".00362"), arb(".01914")
    assert arb(".00311")/arb(".8616")*ratio_errors < ru
    assert arb(".000268")/arb(".01401")*ratio_errors < rd
    # Every matrix is normalized to largest singular value one.
    leakage = arb(".00042")
    assert (1+2*rd**2).sqrt()*arb(".000415") < leakage
    bu, bd = (heavy_line_radius(r, leakage) for r in (ru, rd))
    mix = 1-((arb(".04275")+arb(0,tol)).cos()*(arb(".00377")+arb(0,tol)).cos())**2
    assert bu < arb(".00477") and bd < arb(".02275")
    assert bu+bd < arb(".02751") < arb(".0428") < mix.sqrt()
    assert mix.sqrt() < arb(".043") < 1-(bu+bd)
    # Permit arbitrary additional Yukawa components, measured in the HS norm.
    thick = leakage + arb(".0029")
    radii = sum(heavy_line_radius(r, thick) for r in (ru, rd))
    assert radii < mix.sqrt() and mix.sqrt() < 1-radii
    print("top-bottom projector distance below", (bu+bd).upper().str(9),
          "or above", (1-bu-bd).lower().str(9), "input", mix.sqrt().str(9))
