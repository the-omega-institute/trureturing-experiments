"""Smallest label radius permitting strictly uphill paths from all feasible zeta indices."""
from flint import acb, arb, ctx

FIRST, LAST, TARGET = 10855, 12092, 10875
GUARANTEED = {10871, TARGET} | set(range(10905, 12055))


def scores():
    root_two, log_two = arb(2).sqrt(), arb(2).log()
    error = arb(0, ".000092")
    values = {}
    heights = {}
    for n, z in enumerate(acb.zeta_zeros(FIRST, LAST-FIRST+1), FIRST):
        g = z.imag
        heights[n] = g
        u = root_two*(g*log_two).cos()
        modulus = (3-2*u)*(2*(113+72*u)).sqrt()/13
        denominator = g*g + arb(1)/4
        values[n] = (modulus/denominator, (modulus+error)/denominator)
    return values, heights


def verify():
    ctx.prec = 144
    values, heights = scores()
    peaks = sum(values[n][0] > max(values[n-1][0], values[n+1][0]) for n in range(10906, 12054))
    assert peaks == 213
    assert all(not values[n][0].overlaps(values[n+1][0]) for n in range(10905,12054))
    robust = {n: interval for n, (_, interval) in values.items()}
    # Every possible source has a guaranteed successor with a strict increase,
    # including sources whose own feasibility has not yet been decided.
    margins = []
    for n, current in robust.items():
        if n == TARGET:
            continue
        best_lower = max(robust[j].lower() for j in range(max(FIRST,n-58),min(LAST,n+58)+1)
                         if j in GUARANTEED and j != n)
        assert best_lower > current.upper()
        margins.append(best_lower-current.upper())
    source = 11499
    plateau = max(robust[j].upper() for j in range(source-57,source+58) if j != source)
    assert robust[source].lower() > plateau
    assert robust[11441] > robust[source]
    print("Nearest-neighbor interior strict peaks", peaks)
    print("Sharp uphill radius 58, worst certified edge margin", min(margins).lower().str(9))
    print("Radius 57 obstruction", (robust[source].lower()-plateau).str(9))

    # Chord distance in the actual charge phase w=(g+i/2)/(g-i/2).
    norm = {n: (g*g+arb(1)/4).sqrt() for n, g in heights.items()}

    def chord(n, j):
        return abs(heights[n]-heights[j]) / (norm[n]*norm[j])

    threshold = chord(11499, 11441)
    assert arb("3.92891105968e-7") < threshold < arb("3.92891105969e-7")
    for j, alternative in robust.items():
        if j not in (11499, 11441):
            assert chord(11499,j) > threshold or robust[11499] > alternative

    for n in robust:
        if n in (11499, TARGET):
            continue
        # Within this short list a guaranteed improving successor always exists.
        assert any(j in GUARANTEED and robust[j] > robust[n] and chord(n,j) < threshold
                   for j in range(max(FIRST,n-70), min(LAST,n+70)+1))
    print("Sharp charge-phase chord threshold", threshold.str(16))


if __name__ == "__main__":
    verify()
