"""Exact ordinary even-power abundancy by direct divisor enumeration."""

from dataclasses import dataclass
from fractions import Fraction
from math import isqrt


@dataclass(frozen=True)
class QuadraticRing:
    d: int

    @property
    def a(self):
        return (self.d - 1) // 4 if self.d % 4 == 1 else self.d

    @property
    def b(self):
        return int(self.d % 4 == 1)

    def multiply(self, left, right):
        x, y = left
        u, v = right
        return x * u + self.a * y * v, x * v + y * u + self.b * y * v

    def power(self, z, exponent):
        result = (1, 0)
        for _ in range(exponent):
            result = self.multiply(result, z)
        return result

    def conjugate(self, z):
        x, y = z
        return x + self.b * y, -y

    def norm(self, z):
        x, y = z
        return x * x + self.b * x * y - self.a * y * y

    def in_sector(self, x, y):
        if self.d in {-1, -3}:
            return x > 0 and y >= 0
        return y > 0 or (y == 0 and x > 0)

    def divisor_norms(self, z):
        target = self.norm(z)
        assert target > 0
        # The norm is at least (x² + y²)/2. A divisor norm divides
        # the target norm, so this box includes all divisor representatives.
        limit = isqrt(2 * target)
        norms = []
        for y in range(limit + 1):
            for x in range(-limit, limit + 1):
                if not self.in_sector(x, y):
                    continue
                divisor = (x, y)
                norm = self.norm(divisor)
                if target % norm:
                    continue
                # z/divisor = z * conjugate(divisor) / Norm(divisor).
                u, v = self.multiply(z, self.conjugate(divisor))
                if u % norm == 0 and v % norm == 0:
                    norms.append(norm)
        return sorted(norms)

    def index(self, z, r):
        return Fraction(sum(n ** r for n in self.divisor_norms(z)), self.norm(z) ** r)


def geometric_sum(q, exponent, r):
    return sum(Fraction(1, q ** (r * j)) for j in range(exponent + 1))


def main():
    gaussian = QuadraticRing(-1)
    minus_two = QuadraticRing(-2)
    # Defant, Remark 2.2 and the example after Definition 2.2.
    assert sum(gaussian.divisor_norms((3, 0))) == 10
    assert sum(minus_two.divisor_norms((3, 0))) == 16
    assert gaussian.divisor_norms((9, 3)) == [1, 2, 5, 9, 10, 18, 45, 90]
    assert gaussian.index((9, 3), 1) == 2
    print("Source fidelity: delta_2(3) = 10 at d=-1 and 16 at d=-2; I_2(9+3i) = 2.")

    ring = QuadraticRing(-11)
    pi, theta, w = (0, 1), (1, 1), (-3, 7)
    assert ring.multiply(pi, ring.conjugate(pi)) == (3, 0)
    assert ring.multiply(ring.power(pi, 3), theta) == tuple(-x for x in w)
    assert (ring.norm(pi), ring.norm(theta), ring.norm(w)) == (3, 5, 135)
    assert ring.divisor_norms((3, 0)) == [1, 3, 3, 9]
    assert ring.divisor_norms(w) == [1, 3, 5, 9, 15, 27, 45, 135]
    witness_index = Fraction(16, 9)
    assert ring.index((3, 0), 1) == ring.index(w, 1) == witness_index
    print("Witness at d=-11: norms 9 and 135; delta_2 values 16 and 240; I_2 = 16/9.")

    for e3 in range(4):
        for e5 in range(4):
            lifted = ring.multiply(ring.power(pi, e3), ring.power(theta, e5))
            assert ring.norm(lifted) == 3 ** e3 * 5 ** e5
            for r in range(1, 4):
                expected = geometric_sum(3, e3, r) * geometric_sum(5, e5, r)
                assert ring.index(lifted, r) == expected
    print("Transfer: 48 checks using 0 <= e_3, e_5 <= 3 and 1 <= r <= 3.")

    assert all((a * a - a + 3) % 2 != 0 for a in (0, 1))
    for t in range(6):
        scale = 2 ** t
        first, second = (3 * scale, 0), (-3 * scale, 7 * scale)
        assert (ring.norm(first), ring.norm(second)) == (9 * 4 ** t, 135 * 4 ** t)
        common = witness_index * geometric_sum(4, t, 1)
        assert ring.index(first, 1) == ring.index(second, 1) == common
        print(f"t={t}: norms {ring.norm(first)}, {ring.norm(second)}; I_2 = {common}.")

    print("ALL_CHECKS_OK")


if __name__ == "__main__":
    main()
