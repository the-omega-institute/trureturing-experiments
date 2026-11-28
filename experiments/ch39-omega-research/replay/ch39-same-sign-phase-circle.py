"""Exclude a common phase in the same-sign two-cycle readout for shells seven to nine."""
import argparse
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, arb_mat


ANGLES = ("-2.29931640625", "-0.3143824576973494", "-1.8145357643647222", "0.62744140625",
          "0.484375", "2.5890145395207624", "-2.049072265625", "2.390612743505245",
          "-2.133292208417787", "0.6018105380446326", "0.277587890625", "-0.345947265625")
ACTIVE = (1, 2, 5, 7, 8, 9)


def centered_residual(projected, levels, weights, heat, field):
    low, heavy, top = levels
    newton, vacuum = -arb(1), -arb("1e-122")
    for i, weight in enumerate(weights):
        amplitude = low + (heavy-low)*abs(projected[i,0])**2 + (top-low)*abs(projected[i,1])**2
        lower = max(low.lower(), amplitude.lower())
        upper = min(top.upper(), amplitude.upper())
        assert 0 < lower <= upper
        domain = lower.union(upper)
        center = domain.mid()
        squared = center**2
        k2, k4, slope2, slope4 = (arb(0) for _ in range(4))
        for multiplicity, alpha, mass in field:
            k2 += alpha * heat(2, mass/squared)
            k4 += multiplicity * heat(3, mass/squared)
            exposure = (-mass/domain**2).exp()
            slope2 += alpha * exposure
            slope4 += multiplicity * exposure
        # Sum the signed field derivatives before bounding the input uncertainty.
        shift = domain - center
        newton += weight*(squared*k2/2 + shift*domain*slope2)/arb.pi()
        vacuum -= weight*(squared**2*k4/4 + shift*domain**3*slope4/2)/arb.pi()
    return arb_mat([[newton], [vacuum]])


def verify(root):
    def load(name):
        return run_path(str(root / f"ch39-{name}.py"))

    source = load("fifth-full-positive-joint")
    fifth, columns, levels = source["verify"](root, return_operator="spectral")
    order = load("two-shell-ordered-bond-gates")
    geometry = load("fifth-seventh-shortest-site-gates")["geometry"]
    gravity = load("shortest-common-complex-gate")["gravity_system"]
    station = fifth * columns
    last = order["LAYERS"][-1]
    assert [j for j, _ in order["LAYERS"]] == [0, 1, 2, 2]
    flip = [i for pair in last[1] for i in pair]
    assert len(set(flip)) == len(flip) == 6
    shells = []
    for number in range(2_306_888_100_092, 2_306_888_100_095):
        _, _, _, basis, _, weights, heat, field = geometry(root, fifth, columns, levels, number=number)
        shells.append((basis.conjugate().transpose(), weights, heat, field))

    def system(point, selected, values):
        mats = order["transport"](station, point, parameters=ANGLES, active=ACTIVE)
        # Each disjoint final gate is -G. This is exactly the original circuit
        # followed by a sign on its six final endpoints, for values and derivatives.
        for matrix in mats:
            for i in flip:
                for j in range(matrix.ncols()):
                    matrix[i,j] = -matrix[i,j]
        residuals, derivatives = [], []
        for left, weights, heat, field in selected:
            moved, *slopes = [left*m for m in mats]
            _, jac = gravity(moved, slopes, levels, weights, heat, field, False)
            if values:
                residuals.extend(centered_residual(moved, levels, weights, heat, field).tolist())
            derivatives.extend(jac.tolist())
        return arb_mat(residuals) if values else None, arb_mat(derivatives)


    fixed_parameters = list(ANGLES)
    fixed_parameters[11] = "0"
    first_two = order["transport"](station, [], parameters=fixed_parameters, active=(),
                                   layers=order["LAYERS"][:2])[0]
    entries = order["block"]([arb(x) for x in fixed_parameters[8:]])[0]
    early = order["LAYERS"][2][1]
    late = order["LAYERS"][3][1]
    assert [row[0] for row in early] == [row[0] for row in late]
    assert len(set(x for pair in early for x in pair) | set(x for pair in late for x in pair)) == 9
    a, b, c, d = entries
    coeff = [acb_mat(first_two), acb_mat(first_two.nrows(), first_two.ncols()),
             acb_mat(first_two.nrows(), first_two.ncols())]
    for ((center, early_leaf), (_, late_leaf)) in zip(early, late, strict=True):
        for k in range(first_two.ncols()):
            x, y, z = [first_two[row, k] for row in (center, early_leaf, late_leaf)]
            intermediate = a*x + b*y
            for row in (center, early_leaf, late_leaf):
                coeff[0][row, k] = 0
            coeff[1][center, k], coeff[2][center, k] = b*z, a*intermediate
            coeff[1][early_leaf, k] = c*x + d*y
            coeff[1][late_leaf, k], coeff[2][late_leaf, k] = d*z, c*intermediate

    for angle in ("0", ".37", "-1.29"):
        test = acb(0, arb(angle)).exp()
        parameters = fixed_parameters.copy()
        parameters[11] = angle
        ordinary = order["transport"](station, [], parameters=parameters, active=())[0]
        evaluated = coeff[0] + coeff[1]*test + coeff[2]*(test**2)
        assert all(x.contains(0) for x in (ordinary-evaluated).entries())

    projected = [(tuple(left*matrix for matrix in coeff), weights, heat, field)
                 for left, weights, heat, field in shells[:3]]
    excluded = [0] * 6
    unresolved = []
    deepest = 0
    tested = 0

    def cover(index, depth):
        nonlocal deepest, tested
        tested += 1
        deepest = max(deepest, depth)
        slices = 2**depth
        low, high = (2*arb.pi()*k/slices for k in (index, index+1))
        z = acb(0, low.union(high)).exp()
        for s, (mats, weights, heat, field) in enumerate(projected):
            value = centered_residual(mats[0] + mats[1]*z + mats[2]*(z**2),
                                      levels, weights, heat, field)
            for j in (0, 1):
                if not value[j,0].contains(0):
                    excluded[2*s+j] += 1
                    return
        if depth >= 14:
            unresolved.append((index, depth))
            return
        cover(2*index, depth+1)
        cover(2*index+1, depth+1)

    for index in range(32):
        cover(index, 5)
    print("phase exclusion", "tested", tested, "max depth", deepest,
          "exclusion counts", excluded, "unresolved", len(unresolved), flush=True)
    if unresolved:
        raise SystemExit(f"Unresolved phase intervals: {unresolved[:12]}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("replay_dir", type=Path, help="replay directory at 9e2b8abfe862acb15667d5a0154731a67e491a95")
    args = parser.parse_args()
    if not __debug__:
        parser.error("run without -O because the verification uses assertions")
    verify(args.replay_dir.resolve())
