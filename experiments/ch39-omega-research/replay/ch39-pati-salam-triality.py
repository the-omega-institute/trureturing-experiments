"""Integer charge checks for the split-adjoint Pati–Salam example."""
from itertools import combinations_with_replacement, permutations, product


def verify():
    charges = (1, 1, -2)
    weights = (1, 1, 1, -3)
    # Each row contains X, SU(4) cubic index, and multiplicities of SU(2) doublets.
    families = [(x, 1, 4, 0) for x in charges] + [(-x, -1, 0, 4) for x in charges]
    assert sum(2 * cubic for _, cubic, _, _ in families) == 0
    assert sum(2 * x for x, _, _, _ in families) == 0  # twice the SU(4)^2 X coefficient
    for index in (2, 3):
        assert sum(row[0] * row[index] for row in families) == 0
        assert sum(row[index] for row in families) == 12
    for power in (1, 3):
        assert sum(8 * x**power for x, _, _, _ in families) == 0

    def residual(x, d):
        return (3 * x - 2 * d) % 9

    for x in charges:
        assert [residual(x, d) for d in weights] == [1, 1, 1, 0]
        assert [residual(-x, -d) for d in weights] == [8, 8, 8, 0]
    # A charge-three singlet, the neutral R component, bidoublet and diagonal adjoint.
    assert all(residual(x, d) == 0 for x, d in [(3, 0), (-2, -3), (0, 0)])
    assert {residual(0, d - e) for d, e in product(weights, repeat=2)} == {0, 1, 8}
    assert all((3 * residual(1, d) - center) % 9 == 0
               for d, center in zip(weights, (3, 3, 3, 0)))

    for x, y in product(charges, repeat=2):
        dirac = x - y
        majorana = -x - y - 4
        assert any(dirac + 3 * n == 0 for n in range(-1, 2))
        assert any(majorana + 3 * n == 0 for n in range(3))
    # Integer eigenvalues (up to common denominators) of the colour and R projectors.
    for d, r in product((1, -3), (0, 1)):
        colour, lepton = (d + 3) // 4, (1 - d) // 4
        assert sorted((colour * r, colour * (1-r), lepton * r, lepton * (1-r))) == [0, 0, 0, 1]
    minimal = [xs for xs in combinations_with_replacement(range(-2, 3), 3)
               if sum(xs) == 0 and all(x % 3 == 1 for x in xs)]
    assert minimal == [(-2, 1, 1)]
    # On this monomial spin and weak indices agree. The two epsilon factors
    # on each pair square to one if its indices differ. The colour sign
    # cancels the Grassmann sign of reordering into ascending colour order.
    labels = (0, 1, 0, 1)
    coefficient = sum(labels[a] != labels[b] and labels[c] != labels[d]
                      for a, b, c, d in permutations(range(4)))
    assert coefficient == 16
    print("anomalies, vacuum residues, flavour entries and the QQQL coefficient 16 verified")


if __name__ == "__main__":
    verify()
