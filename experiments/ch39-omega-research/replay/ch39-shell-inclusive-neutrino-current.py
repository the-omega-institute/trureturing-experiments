"""Sharp static flavor-overlap floor for a single local heavy-vector neutrino channel."""
from pathlib import Path
from runpy import run_path
from flint import arb

shell = run_path(str(Path(__file__).with_name('ch39-shell-current-word-floor.py')))
green, points = shell['green'], shell['points']
groups = [{0, 1, 2, 3, 4, 5}, {9, 10, 13, 15, 16, 17}, {11, 12, 18}, {14}]


def main():
    r = green[6, 6].real
    g = green[7, 6]
    h, b = g.real, abs(g.imag)
    center = (r*r - abs(g)**2) / (2 * (r - h))
    radius = r - center
    cutoff = (center - h) / radius
    low = arb('0.9')
    assert 0 < b and low < cutoff < 1
    assert arb('0.99995638902') < cutoff < arb('0.99995638904')

    # Coordinate permutations preserve the Green kernel up to complex conjugation.
    for i in (6, 7, 8):
        other = [j for j in (6, 7, 8) if i != j]
        assert green[i, i].real.overlaps(r)
        assert green[other[0], i].overlaps(green[other[1], i].conjugate())
        for j in range(19):
            if j not in (6, 7, 8):
                assert abs(green[j, i] - center) < radius
                assert green[j, i].real > h

    # Real disk centers give lower bounds on the same real part in every column.
    for group, value in zip(groups, ('.083', '.058', '.071', '.063'), strict=True):
        z = arb(value)
        for i in group:
            rad = max(abs(green[j, i] - z).upper() for j in range(19))
            assert low*z - rad > low*center - radius
            assert cutoff*z - rad > h
    q = arb('.22704').cos() * arb('.00377').cos()
    assert low < q < cutoff
    floor = q*center - radius
    assert arb('0.02075679798') < floor < arb('0.02075679800')
    assert arb('0.00043084466') < floor**2 < arb('0.00043084467')
    mass = (cutoff - q) / (1 + cutoff)
    assert 0 < mass < 1
    print('Transition overlap', cutoff.str(15), 'nominal |Vud|', q.str(15),
          'sharp inclusive amplitude floor', floor.str(16), 'rate floor', (floor**2).str(16))


if __name__ == '__main__':
    main()
