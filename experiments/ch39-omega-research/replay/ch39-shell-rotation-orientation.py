"""Exact orientation test on the nontrivial cyclic-coordinate sector of the shell."""
from pathlib import Path
from runpy import run_path

from sympy import Matrix, eye, zeros

geometry = run_path(str(Path(__file__).with_name("ch39-refined-mass-obstruction.py")))["geometry"]


def orbit_counts(nodes):
    remaining = set(nodes)
    counts = {1: 0, 3: 0}
    while remaining:
        x, y, z = next(iter(remaining))
        orbit = {(x, y, z), (y, z, x), (z, x, y)}
        assert orbit <= remaining
        remaining.difference_update(orbit)
        counts[len(orbit)] += 1
    return counts[1], counts[3]


def verify():
    nodes, _ = geometry()
    assert orbit_counts(nodes[:19]) == (1, 6)
    assert orbit_counts(nodes) == (1, 25)
    # Difference basis e0-e1, e1-e2 on a three-cycle. The reflection swaps e0 and e2.
    gram = Matrix([[2, -1], [-1, 2]])
    rotation = Matrix([[0, -1], [1, -1]])
    flip = Matrix([[0, -1], [-1, 0]])
    assert rotation**2 + rotation + eye(2) == zeros(2)
    assert rotation.T * gram * rotation == gram
    assert flip.T * gram * flip == gram and flip.det() == -1
    commutator = flip * rotation - rotation * flip
    assert commutator.T * gram * commutator == 3 * gram
    print("Original and refined shells: 12 and 50 nontrivial real dimensions.")
    print("Squared orientation defect on either complement: sharp lower bound 3.")


if __name__ == "__main__":
    verify()
