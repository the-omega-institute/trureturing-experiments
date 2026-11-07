"""Finite fixed-point trace for the relational quantum automaton certificate.

The trace lifts the orbit certificate from one final result to an explicit
monotone family S_0 subseteq S_1 subseteq ... and records which finite rules
have passed at each round. It is deliberately a finite proxy. It does not
prove universal successor-kernel or tensor-product closure.
"""
from __future__ import annotations

import json
from pathlib import Path

from verify import EPS, GENERATORS, GATES, apply, canonical_state, norm_squared, probabilities, bloch

RULE_PREREQUISITE = (0, 1, 2, 3, 0, 3)
RULE_NAMES = (
    "R0_identity_norm",
    "R1_binary_split",
    "R2_oriented_ternary",
    "R3_successor_proxy",
    "R4_positive_readout",
    "R5_composition_proxy",
)


def _record(word: str, state: tuple[complex, complex]) -> dict[str, object]:
    p0, p1 = probabilities(state)
    return {
        "word": word or "I",
        "probabilities": [p0, p1],
        "bloch": list(bloch(state)),
        "norm_squared": norm_squared(state),
    }


def _rule_status(round_index: int) -> tuple[list[str], list[str]]:
    passed: list[str] = []
    waiting: list[str] = []
    for name, prerequisite in zip(RULE_NAMES, RULE_PREREQUISITE):
        if round_index >= prerequisite:
            passed.append(name)
        else:
            waiting.append(name)
    return passed, waiting


def closure_trace(max_rounds: int = 12) -> dict[str, object]:
    start = (1.0 + 0.0j, 0.0j)
    states: dict[tuple[float, float, float, float], tuple[str, tuple[complex, complex]]] = {
        canonical_state(start): ("I", start)
    }
    rounds: list[dict[str, object]] = []
    fixed_point = False

    for round_index in range(max_rounds + 1):
        passed, waiting = _rule_status(round_index)
        rounds.append(
            {
                "round": round_index,
                "records": len(states),
                "passed": passed,
                "waiting": waiting,
                "new_records": 0,
            }
        )
        additions: dict[tuple[float, float, float, float], tuple[str, tuple[complex, complex]]] = {}
        for word, state in tuple(states.values()):
            for label in GENERATORS:
                next_state = apply(GATES[label], state)
                key = canonical_state(next_state)
                if key not in states and key not in additions:
                    additions[key] = (word + label if word != "I" else label, next_state)
        rounds[-1]["new_records"] = len(additions)
        if not additions:
            fixed_point = True
            break
        states.update(additions)

    return {
        "seed": ["identity"],
        "rounds": rounds,
        "fixed_point": fixed_point,
        "records": [_record(word, state) for word, state in states.values()],
        "rule_names": list(RULE_NAMES),
    }


def main() -> None:
    output = Path(__file__).with_name("closure_trace.json")
    payload = closure_trace()
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "rounds": len(payload["rounds"]),
        "records": len(payload["records"]),
        "fixed_point": payload["fixed_point"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
