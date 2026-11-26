"""Order the six reported window weights using integer input cells and rational arithmetic."""
from fractions import Fraction

# Lower endpoints for (gamma_n, Re(2^(1/2+i gamma_n))) in units of 2^-20.
CELLS = (
    (10854, 11117846216, -647744),
    (10875, 11136953530, -560739),
    (10843, 11108282848, -698676),
    (10871, 11133461237, -638705),
    (10908, 11165435665, -590558),
    (10918, 11174991557, -544515),
)
SCALE = 2**20


def verify():
    previous_lower = None
    relative_error = Fraction(1, 32768)
    for index, height, real_part in CELLS:
        g, g_max = (Fraction(height + j, SCALE) for j in (0, 1))
        u, u_max = (Fraction(real_part + j, SCALE) for j in (0, 1))
        assert 0 < g and u_max < Fraction(3, 2) and 113 + 72 * u > 0
        low = 2 * (3 - 2*u_max)**2 * (113 + 72*u) / 169
        high = 2 * (3 - 2*u)**2 * (113 + 72*u_max) / 169
        assert low > Fraction(15, 4)**2
        lower = low * (1 - relative_error)**2 / (g_max**2 + Fraction(1, 4))**2
        upper = high * (1 + relative_error)**2 / (g**2 + Fraction(1, 4))**2
        if previous_lower is not None:
            assert previous_lower > upper
        previous_lower = lower
        print(index, "certified")


if __name__ == "__main__":
    verify()
