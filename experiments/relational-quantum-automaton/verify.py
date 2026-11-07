"""Finite relational quantum automaton certificate.

This experiment keeps the smallest complex two-component state space and
checks the closure chain used by the paper:
  unit -> binary split -> three real relational coordinates.
The gates are ordinary one-qubit generators. The experiment records only
finite identities, so it is a regression certificate rather than a proof of
the full Trureturing theory. License: Apache-2.0 under the repository license.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

EPS = 1e-10
State = tuple[complex, complex]
Matrix = tuple[tuple[complex, complex], tuple[complex, complex]]

SQRT2_INV = 1.0 / math.sqrt(2.0)
GATES: dict[str, Matrix] = {
    "I": ((1.0 + 0.0j, 0.0j), (0.0j, 1.0 + 0.0j)),
    "X": ((0.0j, 1.0 + 0.0j), (1.0 + 0.0j, 0.0j)),
    "Z": ((1.0 + 0.0j, 0.0j), (0.0j, -1.0 + 0.0j)),
    "H": (
        (SQRT2_INV + 0.0j, SQRT2_INV + 0.0j),
        (SQRT2_INV + 0.0j, -SQRT2_INV + 0.0j),
    ),
    "S": ((1.0 + 0.0j, 0.0j), (0.0j, 1.0j)),
}
GENERATORS = ("H", "S", "X", "Z")


def matrix_multiply(left: Matrix, right: Matrix) -> Matrix:
    return tuple(
        tuple(
            sum(left[row][k] * right[k][column] for k in range(2))
            for column in range(2)
        )
        for row in range(2)
    )  # type: ignore[return-value]


def apply(matrix: Matrix, state: State) -> State:
    return tuple(
        sum(matrix[row][column] * state[column] for column in range(2))
        for row in range(2)
    )  # type: ignore[return-value]


def norm_squared(state: State) -> float:
    return sum(abs(amplitude) ** 2 for amplitude in state)


def probabilities(state: State) -> tuple[float, float]:
    return (abs(state[0]) ** 2, abs(state[1]) ** 2)


def bloch(state: State) -> tuple[float, float, float]:
    a, b = state
    return (
        2.0 * (a.conjugate() * b).real,
        2.0 * (a.conjugate() * b).imag,
        abs(a) ** 2 - abs(b) ** 2,
    )


def canonical_state(state: State) -> tuple[float, float, float, float]:
    """Quotient a normalized state by its unobservable global phase."""
    pivot = next((amplitude for amplitude in state if abs(amplitude) > EPS), None)
    if pivot is None:
        raise ValueError("zero state")
    phase = pivot / abs(pivot)
    rotated = tuple(amplitude / phase for amplitude in state)
    return tuple(
        round(component, 12)
        for amplitude in rotated
        for component in (amplitude.real, amplitude.imag)
    )


def close_matrix(left: Matrix, right: Matrix, tolerance: float = EPS) -> bool:
    return all(
        abs(left[row][column] - right[row][column]) <= tolerance
        for row in range(2)
        for column in range(2)
    )


def identity_checks() -> dict[str, bool]:
    identity = GATES["I"]
    return {
        "H_squared_is_identity": close_matrix(
            matrix_multiply(GATES["H"], GATES["H"]), identity
        ),
        "S_fourth_power_is_identity": close_matrix(
            matrix_multiply(
                matrix_multiply(GATES["S"], GATES["S"]),
                matrix_multiply(GATES["S"], GATES["S"]),
            ),
            identity,
        ),
        "X_squared_is_identity": close_matrix(
            matrix_multiply(GATES["X"], GATES["X"]), identity
        ),
        "Z_squared_is_identity": close_matrix(
            matrix_multiply(GATES["Z"], GATES["Z"]), identity
        ),
    }


def orbit(max_depth: int = 12) -> list[dict[str, object]]:
    """Enumerate the finite projective orbit of |0> under the generators."""
    start: State = (1.0 + 0.0j, 0.0j)
    queue: list[tuple[str, State]] = [("", start)]
    seen = {canonical_state(start)}
    records: list[dict[str, object]] = []
    while queue:
        word, state = queue.pop(0)
        px, py = probabilities(state)
        x, y, z = bloch(state)
        records.append(
            {
                "word": word or "I",
                "probabilities": [px, py],
                "bloch": [x, y, z],
                "norm_squared": norm_squared(state),
            }
        )
        if len(word) >= max_depth:
            continue
        for label in GENERATORS:
            next_state = apply(GATES[label], state)
            key = canonical_state(next_state)
            if key not in seen:
                seen.add(key)
                queue.append((word + label, next_state))
    return records


def certificate() -> dict[str, object]:
    records = orbit()
    max_norm_error = max(abs(float(row["norm_squared"]) - 1.0) for row in records)
    max_sphere_error = max(
        abs(sum(component * component for component in row["bloch"]) - 1.0)
        for row in records
    )
    max_probability_error = max(
        abs(sum(row["probabilities"]) - 1.0) for row in records
    )
    max_born_error = max(
        max(
            abs(row["probabilities"][0] - (1.0 + row["bloch"][2]) / 2.0),
            abs(row["probabilities"][1] - (1.0 - row["bloch"][2]) / 2.0),
        )
        for row in records
    )
    identities = identity_checks()
    checks = {
        **identities,
        "unit_closure": max_norm_error <= EPS,
        "binary_split_closure": max_probability_error <= EPS,
        "three_coordinate_closure": max_sphere_error <= 5 * EPS,
        "born_readout_identity": max_born_error <= 5 * EPS,
        "projective_orbit_has_six_states": len(records) == 6,
    }
    if not all(checks.values()):
        raise AssertionError(f"certificate failed: {checks}")
    return {
        "model": {
            "state_space": "C^2 modulo global phase",
            "generators": list(GENERATORS),
            "starting_state": "|0> = (1,0)",
            "closure_chain": ["unit", "binary split", "three relational coordinates"],
        },
        "checks": checks,
        "summary": {
            "projective_orbit_size": len(records),
            "max_norm_error": max_norm_error,
            "max_probability_error": max_probability_error,
            "max_bloch_sphere_error": max_sphere_error,
            "max_born_identity_error": max_born_error,
        },
        "states": records,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).with_name("results.json"),
        help="JSON output path (default: results.json beside this script)",
    )
    args = parser.parse_args()
    payload = certificate()
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload["summary"], sort_keys=True))


if __name__ == "__main__":
    main()
