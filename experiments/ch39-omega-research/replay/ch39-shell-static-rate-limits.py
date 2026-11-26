"""Static rates from noisy intensities, pointwise upper limits and zero counts."""
from pathlib import Path
from runpy import run_path
from flint import arb, arb_mat

root = Path(__file__).resolve().parent
tomography = run_path(str(root / 'ch39-shell-arc-rate-tomography.py'))
model = tomography['model']
measure = tomography['measurement_matrix']
support = [1, 2, 3, 4, 5, 6, 9, 12, 17, 23, 33, 49, 77, 130, 233, 410,
           654, 941, 1251, 1570, 1886, 2189, 2469, 2715, 2920, 3073, 3168, 3200]


def signed_rate_design():
    # A fixed dyadic preconditioner keeps the signed dual evaluation accurate.
    matrix = measure(support)
    right = matrix.mid().inv().mid()
    adjusted = matrix * right
    static = measure([0]) * right
    coeff = (static * adjusted.inv()).entries()
    signs = [(-1) ** j for j in range(28)]
    assert all(v * sign > 0 for v, sign in zip(coeff, signs, strict=True))
    factor = sum(abs(v) for v in coeff)
    assert arb('125209.13186765') < factor < arb('125209.13186767')
    dual = adjusted.solve(arb_mat([[s] for s in signs]))
    peak = arb(0)
    for first in range(1, 3201, 64):
        other = [i for i in range(first, min(first + 64, 3201)) if i not in support]
        values = (measure(other) * right) * dual
        for value in values.entries():
            peak = peak.max(abs(value))
    assert arb('0.99999442') < peak < arb('0.99999443')
    print('Mixed rate minimax on indices 1-3200', factor.str(17), 'off-support maximum', peak.str(15))


def coefficients(indices, continuous):
    matrix = arb_mat([model['row'](model['shift'](i), continuous) for i in indices])
    target = arb_mat([model['row'](arb(0), continuous)]).transpose()
    return matrix.transpose().solve(target).entries()


def limits_and_exposure():
    common = [1, 2, 4, 8, 16, 36, 115]
    for continuous, end, lo, hi, cap in [
        (False, 120, '60.71351', '60.71352', '3686.13114551'),
        (True, 115, '61.17696', '61.17697', '3742.62146807'),
    ]:
        indices = common[:-1] + [end]
        # Replay the signed amplitude witness at the prefix and on its entire tail.
        model['verify'](continuous, indices, lo, hi)
        coeff = coefficients(indices, continuous)
        kappa = sum(abs(v) for v in coeff)
        assert arb(cap) < kappa**2 < arb(cap) + arb('1e-8')
        fractions = [abs(v) / kappa for v in coeff]
        print('All-RH null-rate factor', (kappa**2).str(15),
              'zero-event 95 percent factor', (arb(20).log() * kappa**2).str(15))
        print('Exposure fractions', [v.str(8) for v in fractions])
    continuous_coeff = coefficients(common, True)
    total = sum(map(abs, continuous_coeff))
    vertex_coeff = coefficients(common, False)
    shared = sum(c**2 / (abs(b) / total)
                 for c, b in zip(vertex_coeff, continuous_coeff, strict=True))
    assert arb('3688.09472460') < shared < arb('3688.09472462') < total**2
    print('Vertex factor under common exposure', shared.str(16))


if __name__ == '__main__':
    signed_rate_design()
    limits_and_exposure()
