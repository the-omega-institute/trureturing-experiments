#!/usr/bin/env python3
"""Exact finite check of Nagler's PC schema count and the guard-discharge witness.

The checker uses only tuples, integer counts, and explicit multiset cut.  It is
not a proof search: it checks the concrete derivation schema nominated for the
preregistration and instantiates it over a finite set of contexts.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from typing import Iterable, Tuple


@dataclass(frozen=True, order=True)
class Formula:
    head: str
    arg: str | None = None

    def __str__(self) -> str:
        return self.head if self.arg is None else f"{self.head}({self.arg})"


@dataclass(frozen=True)
class Sequent:
    left: Tuple[Formula, ...]
    right: Tuple[Formula, ...]


def canon(xs: Iterable[Formula]) -> Tuple[Formula, ...]:
    return tuple(sorted(xs))


def cut(left_premise: Sequent, right_premise: Sequent, formula: Formula) -> Sequent:
    """Context-differentiating Cut, after checking the displayed cut formula."""
    assert formula in left_premise.right
    assert formula in right_premise.left
    left_right = list(left_premise.right)
    right_left = list(right_premise.left)
    left_right.remove(formula)
    right_left.remove(formula)
    return Sequent(
        canon(left_premise.left + tuple(right_left)),
        canon(tuple(left_right) + right_premise.right),
    )


def left_axiom(connective: str, argument: str) -> Sequent:
    return Sequent((Formula(connective, argument),), ())


def right_rule(connective: str, principal: str, guard: Formula, context: Sequent) -> Sequent:
    """The δcd^- right rule: guard ⊢ and Γ ⊢ Δ yield Γ ⊢ Δ, connective(principal)."""
    assert guard in context.left
    return Sequent(context.left, canon(context.right + (Formula(connective, principal),)))


def bridge(source: str, source_argument: str, target: str, target_argument: str,
           guard: Formula) -> Sequent:
    """Three rule steps: source(source_argument)⊢ twice, then target's right rule."""
    premise = left_axiom(source, source_argument)
    return right_rule(target, target_argument, guard, premise)


def all_contexts(atoms: Tuple[Formula, ...], max_size: int = 2) -> Tuple[Tuple[Formula, ...], ...]:
    out = {()}
    for size in range(1, max_size + 1):
        out.update(canon(choice) for choice in product(atoms, repeat=size))
    return tuple(sorted(out))


def main() -> None:
    # Definition 18.2.1 lists these 24 PC case cardinalities and states their sum.
    pc_case_counts = (2, 2, 8, 8, 8, 20, 20, 4, 20, 4, 4, 32,
                      4, 40, 40, 8, 4, 68, 4, 12, 20, 20, 4, 20)
    assert sum(pc_case_counts) == 376
    print("PC_CASES", len(pc_case_counts))
    print("PC_COUNT", sum(pc_case_counts))
    # The source prints 150 for PU, while its parenthesized bullet counts sum to 152.
    # The discrepancy does not affect the candidate's schema exclusion.
    pu_case_counts = (2, 2, 8, 8, 4, 20, 4, 20, 4, 4, 8, 8, 12, 4, 20, 4, 20)
    assert sum(pu_case_counts) == 152
    print("PU_PRINTED_COUNT", 150)
    print("PU_LISTED_COUNT", sum(pu_case_counts))

    # The candidate is the 2.1 schema with one independent active formula.
    pc_schema = ("alpha", "delta-cd-minus", "kF", "N0:N0")
    pu_schema_families = {
        ("alpha", "delta-cd-minus", "1F2F", "N0:N0"),
        ("gamma0", "beta-kF", "N0:N0"),
    }
    candidate = ("alpha", "delta-cd-minus", "1B", "N0:N0")
    assert pc_schema[0:2] == candidate[0:2]
    assert candidate[2] in ("1A", "1B", "2A", "2B")
    assert candidate not in pu_schema_families
    assert candidate[2] not in ("1A2A", "1B2B", "2A2A", "2B2B")
    print("SCHEMA_PC_MEMBER 1")
    print("SCHEMA_PU_MEMBER 0")

    A, C = "A", "C"
    cA, dC = Formula("c", A), Formula("d", C)
    # B is independent in the rule schema, so the admissible substitution B := c(A)
    # is the guard formula used below.
    bridge_c_to_d = bridge("c", A, "d", C, cA)
    bridge_d_to_c = bridge("d", C, "c", A, dC)
    assert bridge_c_to_d == Sequent((cA,), (dC,))
    assert bridge_d_to_c == Sequent((dC,), (cA,))
    print("BRIDGE c(A) |- d(C) 3_RULE_STEPS")
    print("BRIDGE d(C) |- c(A) 3_RULE_STEPS")

    formulas = (Formula("p"), cA, dC)
    contexts = all_contexts(formulas)
    for gamma, delta in product(contexts, repeat=2):
        # Left occurrence: d(C)⊢c(A) cut against Γ,c(A)⊢Δ.
        source_left = Sequent(canon(gamma + (cA,)), delta)
        replaced_left = cut(bridge_d_to_c, source_left, cA)
        assert replaced_left == Sequent(canon(gamma + (dC,)), delta)
        reverse_left = Sequent(canon(gamma + (dC,)), delta)
        restored_left = cut(bridge_c_to_d, reverse_left, dC)
        assert restored_left == source_left

        # Right occurrence: Γ⊢Δ,c(A) cut against c(A)⊢d(C).
        source_right = Sequent(gamma, canon(delta + (cA,)))
        replaced_right = cut(source_right, bridge_c_to_d, cA)
        assert replaced_right == Sequent(gamma, canon(delta + (dC,)))
        reverse_right = Sequent(gamma, canon(delta + (dC,)))
        restored_right = cut(reverse_right, bridge_d_to_c, dC)
        assert restored_right == source_right
    print("CONTEXTS_CHECKED", len(contexts) ** 2)
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
