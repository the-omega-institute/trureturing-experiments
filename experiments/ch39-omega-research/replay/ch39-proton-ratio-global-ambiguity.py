"""A two-site readout with two different overlap moduli, certified over Q."""

import json
import runpy
from pathlib import Path

from flint import fmpq as Q, fmpq_mat, fmpq_mpoly_ctx


HERE = Path(__file__).resolve().parent
RING = fmpq_mpoly_ctx.get(("U", "V", "rho"))


def readout_pairs():
    source = runpy.run_path(str(HERE / "ch39-proton-two-site-exceptional-orders.py"))
    moments = {}
    for key, poly in source["insertion_moments"]().items():
        # T=t^4=1. All surviving insertion coefficients are nonnegative.
        values = poly.subs({"T": 1}).to_dict()
        assert all(c > 0 for c in values.values())
        moments[key] = RING.from_dict({m[:3]: Q(c) for m, c in values.items()})
    z, a, b = (moments[key] for key in ((0, 0), (1, 0), (0, 1)))
    return ((z*moments[1, 1], a*b), (z*moments[2, 0], a*a), (z*moments[0, 2], b*b))


def jacobian_at(equations, point):
    return fmpq_mat([[poly.derivative(j)(*point) for j in range(3)] for poly in equations])


def coordination_signs():
    data = json.loads((HERE / "ch39-proton-overlap-jacobian-coordination-law.json").read_text())
    terms = [(int(key.split(",")[0]), int(key.split(",")[1]), int(value))
             for key, value in data["P"].items()]
    w = Q(1, 2)
    signs = [sum(c*w**a*r**b for a, b, c in terms) for r in (Q(1, 2), Q(3, 4))]
    assert Q(data["C"]) > 0 and signs[0] > 0 > signs[1]
    print("P(1/2,rho) has opposite signs at rho=1/2 and rho=3/4")


def main():
    coordination_signs()
    pairs = readout_pairs()
    first = (Q(1, 4), Q(1, 8), Q(1, 2))
    target = [num(*first)/den(*first) for num, den in pairs]
    equations = [num - value*den for (num, den), value in zip(pairs, target)]
    assert all(poly(*first) == 0 for poly in equations)
    assert jacobian_at(equations, first).det() != 0
    print("Rational target (AB, AA, BB):", *target)

    center = tuple(Q(i, 10**16) for i in (1476388875849972, 560800415954564, 7716008713962137))
    radius = Q(1, 10**11)
    upper = [x + radius for x in center]
    assert 0 < center[1]-radius and upper[1] < center[0]-radius
    assert Q(77, 100) < center[2]-radius < upper[2] < Q(78, 100)

    inverse = jacobian_at(equations, center).inv()
    residual = inverse * fmpq_mat([[poly(*center)] for poly in equations])
    displacement = max(abs(residual[i, 0]) for i in range(3))

    # For positive arguments every Hessian entry of F_i is bounded on
    # the cube by the same derivative of the absolute-coefficient polynomial
    # evaluated at the upper corner. Sum over both derivative indices.
    hessian_bounds = []
    for poly in equations:
        majorant = RING.from_dict({m: abs(c) for m, c in poly.to_dict().items()})
        hessian_bounds.append(sum(majorant.derivative(j).derivative(k)(*upper)
                                  for j in range(3) for k in range(3)))
    contraction = radius*max(sum(abs(inverse[i, j])*hessian_bounds[j] for j in range(3))
                             for i in range(3))
    assert displacement < radius/10 and contraction < Q(1, 10)
    print("Infinity-norm Newton displacement < radius/10; contraction bound < 1/10")
    print("Second preimage lies within 1e-11 of", *center)
    print("Both preimages are regular, with 0<V/U<1 and disjoint overlap intervals.")


if __name__ == "__main__":
    main()
