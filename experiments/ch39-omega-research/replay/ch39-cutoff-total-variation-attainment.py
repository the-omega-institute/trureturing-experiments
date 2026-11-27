"""Certify the closed-domain total-variation optimum for the three visible spectra.

Run from a Git checkout with Python 3.12+ and python-flint 0.8.0.
The scripts and data used as inputs are extracted from the pinned Git tree.

Krawczyk inclusion proves simultaneous moment matching and dual contact.
The derivative bounds give exactly four dual contacts with signs −,+,−,+.
Contact-column independence forces every optimal measure to contain the
negative atom at zero. To approach the optimum with positive cutoffs, move
that atom to a small positive cutoff and use the positive basis to correct
the moment residual. The correction tends to zero.

For a small positive cutoff floor, the invertible primal–dual Jacobian
continues the solution with its leftmost atom at the floor. Its dual
coefficients stay in the certified box for sufficiently small floors.
The right derivative of the optimal cost is minus the leftmost weight
times the derivative of the dual function at zero. The printed enclosure
uses the derivative of the massless term in the Newton heat trace.

Pass --floor-window to certify the unique four-atom optimum uniformly for
0 <= delta <= 10^-14, with its first atom at delta, for all three spectra.
This mode also certifies the strict convexity of the optimal cost, a quadratic remainder, and
full rank for the contact and tangent columns of the signed Jacobian.
"""

from argparse import ArgumentParser
from pathlib import Path
from math import floor
from runpy import run_path
import subprocess
import sys
import tarfile
import tempfile
from flint import acb, arb, arb_mat, ctx

INPUT_COMMIT = "0382ed8bb2cd901b85f322063845a616acca1b88"
REPLAY_PATH = "experiments/ch39-omega-research/replay"
ENTRY = "ch39-stationary-cutoff-total-variation.py"
GRID_TIMES = {
    11661: ("-12.394768334536108", "-1.0978622593626126"),
    11934: ("-12.470609962039484", "-1.2424734994092408"),
    12207: ("-12.391938402374743", "-1.0937896492338242"),
}
W_SCALE = list(map(arb, ("160000", "160000", "2", "0.25")))
BASE_WEIGHTS = list(map(arb, ("-1", "1", "-1", "1")))


def interval_jets(heat, rows, t):
    if t is None:
        return [arb(1), *[arb(0)] * 4], [arb(0)] * 5, [arb(0)] * 5
    lo, hi = arb(t.lower()), arb(t.upper())
    ph, ps, ell, p, q = heat(rows, lo, hi)
    denom = 2 * arb.pi()
    return (
        [arb(1), ell, ph[0], ps[0], p],
        [arb(0), ph[1] / denom, ph[0] + ph[1], 2 * ps[0] + ps[1], p + q],
        [
            arb(0),
            ph[2] / denom,
            ph[0] + ph[1] + ph[2],
            4 * ps[0] + 3 * ps[1] + ps[2],
            p + q + ps[2] / denom,
        ],
    )


def dot(xs, ys):
    return sum(x * y for x, y in zip(xs, ys))


def identity(n):
    return arb_mat([[arb(int(i == j)) for j in range(n)] for i in range(n)])


def endpoint_column(rows, floor_upper):
    """Enclose the five moments of one cutoff over [0, floor_upper]."""
    x = arb(0).union(floor_upper)
    ell = phi = psi = pressure = arb(0)
    for g, a, z in rows:
        if z.is_zero():
            phi += a * x
            psi += g * x**2
            continue
        scaled_mass = z / floor_upper
        assert scaled_mass >= 1
        # For u >= 1 and n >= 1, E_n(u) <= exp(-u)/u <= exp(-u).
        factor = arb(0).union((-scaled_mass).exp())
        ell += z * a * factor
        phi += a * x * factor
        psi += g * x**2 * factor
        pressure += z * g * x * factor
    denom = 2 * arb.pi()
    return [arb(1), ell / denom, phi, psi, pressure / denom]


def endpoint_derivatives(rows, floor_upper):
    """Enclose the first two cutoff derivatives, including their limits at zero."""
    x = arb(0).union(floor_upper)
    first = [arb(0) for _ in range(5)]
    second = first.copy()
    denom = 2 * arb.pi()
    for g, a, z in rows:
        if z.is_zero():
            first[2] += a
            first[3] += 2 * g * x
            second[3] += 2 * g
            continue
        assert z > 0
        lower = (z / floor_upper).lower()
        # The maximum of u^k exp(-u) on [lower, infinity) occurs at max(lower, k).
        tails = []
        for k in range(4):
            point = max(lower, arb(k))
            tails.append(arb(0).union(point**k * (-point).exp()))
        e, ue, u2e, u3e = tails
        first[1] += a * ue / denom
        first[2] += a * (e + ue)
        first[3] += g * x * (2 * e + ue)
        first[4] += g * z * e / denom
        second[1] += a * (u3e - u2e) / (z * denom)
        second[2] += a * u3e / z
        second[3] += g * (2 * e + 2 * ue + u2e)
        second[4] += g * u2e / denom
    return first, second


def solve(certificate, rows, n, *, floor_upper=None, radius="1e-30"):
    heat = certificate["heat"]
    yscale = [abs(arb(v)) for v in certificate["DUAL"][n]]
    yseed = [arb(v) / s for v, s in zip(certificate["DUAL"][n], yscale)]
    times = [arb(v) for v in GRID_TIMES[n]]
    scales = list(map(arb, ("320000", "2", "0.4", "100", "1e-5")))
    desired = [arb(1), arb(-1), arb(0), arb(0), 4 * arb("1e-122")]
    signs = BASE_WEIGHTS

    def system(u, first_column=None):
        dual = [v * s for v, s in zip(u[:5], yscale)]
        weights = [v * s for v, s in zip(u[5:9], W_SCALE)]
        vals, ds, dds = (
            list(z)
            for z in zip(
                *(interval_jets(heat, rows, t) for t in [None, *u[9:], arb(0)])
            )
        )
        if first_column is not None:
            vals[0] = first_column
        f = [
            (dot([row[i] for row in vals], weights) - desired[i]) / scales[i]
            for i in range(5)
        ]
        f += [dot(v, dual) - s for v, s in zip(vals, signs)]
        f += [dot(row, dual) for row in ds[1:3]]
        jac = [[arb(0) for _ in range(11)] for _ in range(11)]
        for i in range(5):
            for j in range(4):
                jac[i][5 + j] = vals[j][i] * W_SCALE[j] / scales[i]
            for j in range(2):
                jac[i][9 + j] = ds[1 + j][i] * weights[1 + j] / scales[i]
        for i in range(4):
            for j in range(5):
                jac[5 + i][j] = vals[i][j] * yscale[j]
        for i in range(2):
            jac[6 + i][9 + i] = dot(ds[1 + i], dual)
            jac[9 + i][9 + i] = dot(dds[1 + i], dual)
            for j in range(5):
                jac[9 + i][j] = ds[1 + i][j] * yscale[j]
        return arb_mat(11, 1, f), arb_mat(jac)

    u = [*yseed, *BASE_WEIGHTS, *times]
    for _ in range(8):
        f, J = system(u)
        correction = J.solve(f)
        u = [(value - correction[i, 0]).mid() for i, value in enumerate(u)]
    _, J = system(u)
    r = arb(radius)
    rho = arb(0, r)
    box = [value + rho for value in u]
    endpoint = None if floor_upper is None else endpoint_column(rows, floor_upper)
    f, _ = system(u, endpoint)
    _, B = system(box, endpoint)
    inv = J.inv()
    inverse = arb_mat([[inv[i, j].mid() for j in range(11)] for i in range(11)])
    contraction = identity(11) - inverse * B
    row_bounds = [
        sum(abs(contraction[i, j]) for j in range(11)).upper() for i in range(11)
    ]
    assert all(bound < 1 for bound in row_bounds)
    contraction_norm = max(row_bounds)
    err = -inverse * f + contraction * arb_mat(11, 1, [rho] * 11)
    if not all(abs(err[i, 0]) < r for i in range(11)):
        print(
            "Krawczyk step",
            [err[i, 0].str(4) for i in range(11)],
            "row norm",
            contraction_norm,
        )
        print("center system", [f[i, 0].str(4) for i in range(11)])
        raise ArithmeticError(n)
    y = [v * s for v, s in zip(box[:5], yscale)]
    w = [v * s for v, s in zip(box[5:9], W_SCALE)]
    print(
        n,
        "Krawczyk radius",
        r.str(4),
        "Lipschitz bound",
        str(contraction_norm),
        "log nodes",
        *[v.str(31) for v in box[9:]],
        "TV enclosure",
        sum(map(abs, w)).str(25),
    )
    assert all(wi * si > 0 for wi, si in zip(w, signs))
    if floor_upper is not None:
        first, second = endpoint_derivatives(rows, floor_upper)
        # Only the five moment equations and the left contact depend directly on the floor.
        parameter_derivative = arb_mat(
            11, 1,
            [w[0] * value / scale for value, scale in zip(first, scales)]
            + [dot(y, first)] + [arb(0)] * 5,
        )
        response = inverse * parameter_derivative
        response_norm = max(abs(response[i, 0]).upper() for i in range(11))
        derivative_norm = response_norm / (1 - contraction_norm)
        # Intersect each enclosure of the derivative with its fixed-point image.
        derivatives = arb_mat(11, 1, [arb(0, derivative_norm)] * 11)
        for _ in range(8):
            image = -response + contraction * derivatives
            derivatives = arb_mat(11, 1, [
                derivatives[i, 0].intersection(image[i, 0]) for i in range(11)
            ])
        dual_derivatives = [derivatives[i, 0] * scale for i, scale in enumerate(yscale)]
        weight_derivative = derivatives[5, 0] * W_SCALE[0]
        second_derivative = -(
            weight_derivative * dot(y, first)
            + w[0] * (dot(dual_derivatives, first) + dot(y, second))
        )
        assert arb("1.13e14") < second_derivative and second_derivative < arb("1.19e14")
        print(n, "optimal cost second derivative", second_derivative.str(15))
        remainder = abs(second_derivative).upper() / 2
        print(n, "quadratic remainder coefficient upper bound", remainder.str(15))
        assert remainder < arb("6e13")
    return y, w, box[9:]


def hermitian_inertia(matrix, shift):
    """Count negative eigenvalues by Hermitian Schur complements."""
    rows = matrix.tolist()
    negatives = 0
    for i in range(len(rows)):
        pivot = rows[i][i].real - shift
        if pivot.contains(0):
            return None
        negatives += int(pivot < 0)
        for j in range(i + 1, len(rows)):
            ratio = rows[j][i] / pivot
            for k in range(i + 1, j + 1):
                rows[j][k] -= ratio * rows[k][i].conjugate()
    return negatives


def hermitian_levels(matrix):
    bit_width = arb(2) ** -156
    low, high = arb(-1), arb(64)
    assert (
        hermitian_inertia(matrix, low) == 0
        and hermitian_inertia(matrix, high) == matrix.nrows()
    )

    def brackets(lo, hi, first, last):
        if first == last:
            return []
        if hi - lo < bit_width:
            return [lo.union(hi)] * (last - first)
        for alpha in (arb(1) / 2, arb(3) / 8, arb(5) / 8):
            mid = (lo + alpha * (hi - lo)).mid()
            count = hermitian_inertia(matrix, mid)
            if count is not None:
                assert first <= count <= last
                return brackets(lo, mid, first, count) + brackets(mid, hi, count, last)
        raise ArithmeticError("all shifted pivots inconclusive")

    return brackets(low, high, 0, matrix.nrows())


def exact_spectra(certificate):
    previous = certificate["previous"]
    old = previous["spectra"]()
    minimal = previous["minimal"]
    supplier = minimal["chiral"]["source"]
    model = supplier["mean"]["inputs"]()[0]
    old_shell = model["shell"]
    smooth = old_shell["smooth"]
    soft = smooth["soft"]
    matter = smooth["matter"]
    assert soft["geo"]["sites"] == supplier["geometry"]["sites"]
    ctx.prec = 320
    assert old.keys() == GRID_TIMES.keys()
    all_rows = {}
    for n, coarse in old.items():
        height = acb.zeta_zero(n).imag
        base, neutral, _ = soft["vector_pencil"](height)
        removal = (
            arb(23)
            / 25
            * arb(2) ** -20
            / (height**2 + soft["point"]["study"]["offset"])
        )
        spoke = (1 - removal) * base
        modes = hermitian_levels(spoke + neutral)
        assert modes[0].contains(0) and all(v > 0 for v in modes[1:])
        # The neutral graph Laplacian has constant vector in its kernel.
        neutral_levels = modes[1:]
        assert len(neutral_levels) == 18
        center = sum(neutral_levels) / 18
        layers = {
            0: minimal["chiral"]["ordered"](
                [v * abs(v - center) for v in neutral_levels]
            )
        }
        for q in model["charges"]:
            eigs = (spoke + matter["nonspokes"](q, height)).eig(multiple=True)
            assert len(eigs) == 19 and all(
                v.real > 0 and v.imag.contains(0) for v in eigs
            )
            layers[q] = minimal["chiral"]["ordered"](
                [v.real * abs(v.real - center) for v in eigs]
            )
        exact = minimal["fields"](layers)
        assert len(exact) == len(coarse)
        mismatches = [
            (i, g, G, z.str(24), Z.str(24), (z - Z).str(9))
            for i, ((g, a, z), (G, A, Z)) in enumerate(zip(exact, coarse))
            if not (g == G and a.overlaps(A) and Z.contains(z))
        ]
        if mismatches:
            print(n, "mismatch count", len(mismatches), "first", mismatches[:7])
            raise ArithmeticError("full spectrum outside pinned enclosure")
        assert all(z.is_zero() or z > 0 for _, _, z in exact)
        all_rows[n] = exact
        print(
            n, "maximum squared-mass radius", max(z.rad() for _, _, z in exact).str(6)
        )
    return all_rows


def derivative_tail_positive(y, rows):
    X = (-arb(40)).exp()
    x = arb(0).union(X)
    error = arb(0)
    leading = arb(0)
    for g, a, z in rows:
        if z.is_zero():
            leading += y[2] * a + 2 * x * y[3] * g
            continue
        lo = arb((z / X).lower())
        assert lo > 2
        tail = lo * lo * (-lo).exp()
        # For u >= lo > 2, the two normalized polynomial factors
        # (1+u) and (2+u) are bounded by u^2. Also E_n(u) <= e^(-u).
        error += tail * (
            abs(y[1] * a) / (2 * arb.pi())
            + abs(y[2] * a)
            + X * abs(y[3] * g)
            + X * abs(y[4] * g) / arb.pi()
        )
    assert leading - error > 0
    return (leading - error).lower()


def third_derivative(heat, rows, y, box):
    lo, hi = arb(box.lower()), arb(box.upper())
    phi, psi, _, p, q = heat(rows, lo, hi)
    x = lo.exp().union(hi.exp())
    phi3 = psi3 = arb(0)
    for g, a, z in rows:
        if z.is_zero():
            continue
        u_lo, u_hi = (z / hi.exp()).lower(), (z / lo.exp()).upper()
        e = (-u_hi).exp().union((-u_lo).exp())
        psi3 += g * z**3 * e / x
        phi3 += a * z**3 * e / x**2
    denom = 2 * arb.pi()
    return (
        y[1] * phi3 / denom
        + y[2] * (phi[0] + phi[1] + phi3)
        + y[3] * (8 * psi[0] + 7 * psi[1] + 3 * psi[2] + psi3)
        + y[4] * (p + q + (psi[2] + psi3) / denom)
    )


def certify_dual(certificate, rows, y, times, *, floor_upper=0):
    heat = certificate["heat"]
    a, b = [floor(float(t) * 1024) for t in times]
    l1, u1, l2, u2 = [arb(k) / 1024 for k in (a - 8, a + 9, b - 8, b + 9)]
    assert -40 < l1 and floor_upper < l1.exp()
    assert times[0] > l1 and times[0] < u1 and times[1] > l2 and times[1] < u2
    assert u1 < l2 and u2 < 0
    count = [0, 0]

    def covering(lo, hi, sign, curvature, depth=0):
        count[0] += 1
        box = lo.union(hi)
        radius = (hi - lo) / 2
        _, d1, d2 = interval_jets(heat, rows, (hi + lo) / 2)
        d2y, d3y = dot(d2, y), third_derivative(heat, rows, y, box)
        if curvature:
            accepted = sign * d2y > radius * abs(d3y)
        else:
            accepted = sign * dot(d1, y) > radius * (abs(d2y) + radius * abs(d3y) / 2)
        if accepted:
            return
        assert depth < 25, (lo.str(12), hi.str(12), sign, curvature)
        count[1] += 1
        mid = ((hi + lo) / 2).mid()
        covering(lo, mid, sign, curvature, depth + 1)
        covering(mid, hi, sign, curvature, depth + 1)

    for lo, hi, positive, curvature in [
        (arb(-40), l1, 1, False),
        (l1, u1, -1, True),
        (u1, l2, -1, False),
        (l2, u2, 1, True),
        (u2, arb(0), 1, False),
    ]:
        while lo < hi:
            endpoint = lo + arb(1) / 16
            if endpoint > hi:
                endpoint = hi
            covering(lo, endpoint, positive, curvature)
            lo = endpoint
    return count


def certify_rank(certificate, rows, times, *, endpoint=None):
    heat = certificate["heat"]
    atoms = [interval_jets(heat, rows, t)[0] for t in [None, *times, arb(0)]]
    if endpoint is not None:
        atoms[0] = endpoint
    minors = [arb_mat([r[:i] + r[i + 1 :] for r in atoms]).det() for i in range(5)]
    assert any(not v.contains(0) for v in minors)
    basis = [interval_jets(heat, rows, arb(t))[0] for t in (-32, -23, -13, -2, 0)]
    det = arb_mat(basis).det()
    assert not det.contains(0)
    return max(abs(v).lower() for v in minors), det


def main():
    if not __debug__ or sys.version_info < (3, 12):
        sys.exit("Use Python 3.12 or newer with assertions enabled.")

    parser = ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--floor-window", action="store_true")
    options = parser.parse_args()

    repo = Path(__file__).resolve().parents[3]
    with tempfile.TemporaryDirectory(prefix="cutoff-attainment.") as directory:
        scratch = Path(directory)
        archive = scratch / "inputs.tar"
        subprocess.run(
            [
                "git",
                "archive",
                "--format=tar",
                "-o",
                str(archive),
                INPUT_COMMIT,
                REPLAY_PATH,
            ],
            cwd=repo,
            check=True,
        )
        with tarfile.open(archive) as inputs:
            inputs.extractall(scratch, filter="data")
        certificate = run_path(str(scratch / REPLAY_PATH / ENTRY))
        data = exact_spectra(certificate)
        ctx.prec = 320
        for n, rows in data.items():
            if options.floor_window:
                floor_upper = arb("1e-14")
                endpoint = endpoint_column(rows, floor_upper)
                y, _, times = solve(certificate, rows, n, floor_upper=floor_upper, radius="1e-9")
                tail = derivative_tail_positive(y, rows)
                intervals = certify_dual(certificate, rows, y, times, floor_upper=floor_upper)
                rank, _ = certify_rank(certificate, rows, times, endpoint=endpoint)
                jets = [interval_jets(certificate["heat"], rows, t) for t in times]
                assert dot(jets[0][2], y) < 0 and dot(jets[1][2], y) > 0
                final = interval_jets(certificate["heat"], rows, arb(0))[0]
                contacts = [endpoint, jets[0][0], jets[1][0], final]
                tangent_minors = [arb_mat(contacts + [jet[1]]).det() for jet in jets]
                assert any(not det.contains(0) for det in tangent_minors)
                print(n, "contact-plus-tangent determinants",
                      *(det.str(9) for det in tangent_minors))
                print(
                    n, "uniform floor window [0, 1e-14]", intervals,
                    "tail lower bound", tail.str(6), "contact minor", rank.str(6),
                )
                continue
            y, weights, times = solve(certificate, rows, n)
            tail = derivative_tail_positive(y, rows)
            intervals = certify_dual(certificate, rows, y, times)
            rank, basis_det = certify_rank(certificate, rows, times)
            massless_newton_weight = sum(a for _, a, z in rows if z.is_zero())
            floor_slope = -weights[0] * y[2] * massless_newton_weight
            assert floor_slope > 0
            print(
                n,
                "strict dual shape",
                intervals,
                "tail lower bound for F derivative divided by x",
                tail.str(6),
                "contact minor",
                rank.str(6),
                "positive basis determinant",
                basis_det.str(6),
                "optimal cost right derivative at zero floor",
                floor_slope.str(25),
            )


if __name__ == "__main__":
    main()
