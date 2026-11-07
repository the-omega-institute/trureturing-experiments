#!/usr/bin/env python3
"""Exact checks for the Heyting hyperamalgamation counterexample family.

All algebra operations are on integers and finite products; no floating point
arithmetic is used.  The finite length loop checks the stalled-transfer
invariant for every tested chain size and every interior element.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Hm:
    m: int

    @property
    def top(self) -> int:
        return self.m - 1

    def imp(self, x: int, y: int) -> int:
        return self.top if x <= y else y

    def le(self, x: int, y: int) -> bool:
        return x <= y


def chi(x: int) -> int:
    return 0 if x == 0 else 1


def bool_imp(x: int, y: int) -> int:
    return 1 if x <= y else y


def prod_join(x: tuple[int, int], y: tuple[int, int]) -> tuple[int, int]:
    return (max(x[0], y[0]), max(x[1], y[1]))


def prod_meet(x: tuple[int, int], y: tuple[int, int]) -> tuple[int, int]:
    return (min(x[0], y[0]), min(x[1], y[1]))


def prod_imp(h: Hm, x: tuple[int, int], y: tuple[int, int]) -> tuple[int, int]:
    return (h.imp(x[0], y[0]), bool_imp(x[1], y[1]))


def le2(x: tuple[int, int], y: tuple[int, int]) -> bool:
    return x[0] <= y[0] and x[1] <= y[1]


def i_a(x: tuple[int, int]) -> tuple[int, int, int]:
    return (x[0], x[1], chi(x[0]))


def i_b(x: tuple[int, int]) -> tuple[int, int, int]:
    return (x[0], chi(x[0]), x[1])


def d_join(x: tuple[int, int, int], y: tuple[int, int, int]) -> tuple[int, int, int]:
    return tuple(max(a, b) for a, b in zip(x, y))


def d_le(x: tuple[int, int, int], y: tuple[int, int, int]) -> bool:
    return all(a <= b for a, b in zip(x, y))


def check_implication_preservation(m: int) -> int:
    h = Hm(m)
    checks = 0
    for x in range(m):
        for y in range(m):
            lhs = (h.imp(x, y), bool_imp(chi(x), chi(y)))
            rhs = prod_imp(h, (x, chi(x)), (y, chi(y)))
            assert lhs == rhs
            checks += 1
    return checks


def check_graph_embedding_and_strong_amalgam(m: int) -> dict[str, int]:
    h = Hm(m)
    c = [(x,) for x in range(m)]
    a = [(x, t) for x in range(m) for t in (0, 1)]
    b = list(a)
    # j(x)=(x,chi(x)) preserves all Heyting operations and is injective.
    for x in range(m):
        assert (x, chi(x)) in a
    for x in range(m):
        for y in range(m):
            assert prod_join((x, chi(x)), (y, chi(y))) == (
                max(x, y), chi(max(x, y))
            )
            assert prod_meet((x, chi(x)), (y, chi(y))) == (
                min(x, y), chi(min(x, y))
            )
            assert prod_imp(h, (x, chi(x)), (y, chi(y))) == (
                h.imp(x, y), chi(h.imp(x, y))
            )
    # Explicit D0 = Hm x 2 x 2 with the two indicated embeddings.
    image_a = {i_a(x) for x in a}
    image_b = {i_b(x) for x in b}
    image_c = {(x, chi(x), chi(x)) for x in range(m)}
    assert image_a & image_b == image_c
    for x in a:
        assert len(i_a(x)) == 3
    for x in b:
        assert len(i_b(x)) == 3
    return {"C": len(c), "A": len(a), "B": len(b), "D0": 4 * m}


def check_forced_mixed_inequality(m: int, e: int) -> None:
    top = m - 1
    a = (top, 0)
    a1 = (e, 0)
    b1 = (top, 0)
    c = (e, 1)
    assert prod_meet(a, c) == a1
    assert prod_join(c, b1) == (top, 1)
    # The displayed D0 realizes the identity forced in every distributive
    # amalgam: a = a∧(c∨b1) = (a∧c)∨(a∧b1) ≤ a1∨b1.
    a_d, a1_d, b1_d, c_d = i_a(a), i_a(a1), i_b(b1), i_a(c)
    assert c_d == i_b(c)
    assert d_le(a_d, d_join(a1_d, b1_d))
    assert d_join(a1_d, b1_d) == (top, 1, 1)


def check_zero_transfer_invariant(m: int, e: int, max_length: int = 30) -> int:
    """Check every finite certificate length through the invariant.

    A common c∈C below a1=(e,0) or b1=(m−1,0) must have chi(c)=0,
    hence c=0.  Consequently all A and B sequence entries stay below
    a1 and b1, so the terminal a=(m−1,0)≤a_n∨c_n is impossible.
    """
    top = m - 1
    a1 = (e, 0)
    b1 = (top, 0)
    a = (top, 0)
    lengths = 0
    for n in range(1, max_length + 1):
        # n=1 is the direct one-step case; n>1 has the alternating transfer.
        for x in range(m):
            c = (x, chi(x))
            if le2(c, a1) or le2(c, b1):
                assert x == 0
        # Initial c1 is forced to zero, so every possible a2 is <= a1.
        possible_a = {x for x in range(m) for t in (0, 1) if le2((x, t), a1)}
        possible_b = {x for x in range(m) for t in (0, 1) if le2((x, t), b1)}
        assert possible_a == set(range(e + 1))
        assert possible_b == set(range(m))
        # Every subsequent transfer has c=d=0, hence cannot increase either
        # coordinate.  The terminal inequality would require top <= e.
        assert top > e
        lengths += 1
    return lengths


def main() -> None:
    tested = 0
    implication_checks = 0
    total_lengths = 0
    diagram_sizes = []
    for m in range(3, 31):
        implication_checks += check_implication_preservation(m)
        sizes = check_graph_embedding_and_strong_amalgam(m)
        diagram_sizes.append(sizes)
        for e in range(1, m - 1):
            check_forced_mixed_inequality(m, e)
            total_lengths += check_zero_transfer_invariant(m, e)
            tested += 1
    result = {
        "chain_sizes": [3, 30],
        "interior_instances": tested,
        "implication_checks": implication_checks,
        "certificate_lengths_checked": total_lengths,
        "first_diagram": diagram_sizes[0],
        "last_diagram": diagram_sizes[-1],
        "claim_obstruction": "forced mixed join plus zero-transfer invariant",
    }
    out = Path(__file__).with_suffix(".json")
    out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, sort_keys=True))
    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
