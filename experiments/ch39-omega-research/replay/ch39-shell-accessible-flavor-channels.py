"""Sharp joint static bounds for two neutrino flavor rows and the charged-lepton row."""
from pathlib import Path
from runpy import run_path
from flint import arb

shell = run_path(str(Path(__file__).with_name('ch39-shell-current-word-floor.py')))


def main():
    green = shell['green']
    minimum = green[7, 6].real
    orbit = (6, 7, 8)
    for i in range(19):
        for j in range(19):
            if i in orbit and j in orbit and i != j:
                assert green[i, j].real.overlaps(minimum)
            else:
                assert green[i, j].real > minimum
    assert green[7, 6].overlaps(green[8, 6].conjugate())

    t12, t13 = arb('.22704'), arb('.00377')
    a, b, c = t12.cos()*t13.cos(), t12.sin()*t13.cos(), t13.sin()
    frame_min = c**2 / (1 + (1 - c**2).sqrt())
    max_row_min = c / (2*c**2 + (1 + a + b)**2).sqrt()
    total = minimum**2 * frame_min
    largest = (minimum * max_row_min)**2
    assert arb('3.4654713929e-9') < total < arb('3.4654713930e-9')
    assert arb('1.4327510507e-9') < largest < arb('1.4327510509e-9')
    print('Three coefficient rates, sharp sum', total.str(16),
          'sharp largest', largest.str(16))


if __name__ == '__main__':
    main()
