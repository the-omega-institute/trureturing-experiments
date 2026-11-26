"""Inner and outer height ranges when the gauge/Planck matching scale is allowed to vary."""
from fractions import Fraction as F
from pathlib import Path
from runpy import run_path
from flint import arb

ROOT = Path(__file__).resolve().parent
fixed = run_path(str(ROOT/'ch39-high-zero-planck-window.py'))
model, q, ratios, offset = [fixed[k] for k in ('model','q','V','OFFSET')]
STEPS = 256


def track(box):
    data = []
    def record(j, inverse, deriv):
        if j >= 23*STEPS:
            data.append([x/q(v) for x, v in zip(inverse, ratios)])
    model['running']['integrate'](box, 33*STEPS, record, subdivisions=STEPS)
    assert len(data) == 10*STEPS
    return data


def verify():
    centers = model['GAUGE']
    bounds = [(q(c*F(999,1000)), q(c*F(1001,1000))) for c in centers]
    box = [a.union(b) for a, b in bounds]
    all_values = track(box)
    low, high = [], []
    for i in range(3):
        for side, dest in enumerate((low, high)):
            face = box.copy()
            face[i] = bounds[i][side]
            dest.append([r[i] for r in track(face)])
    coeff = 4*(fixed['geo']['BOOK_PLANCK']/1000)**2*(3*arb(5).sqrt()/2)/(8*arb.pi())
    possible, sufficient = [], []
    for k, targets in enumerate(all_values):
        start = q(F(23*STEPS+k, STEPS))
        stop = q(F(23*STEPS+k+1, STEPS))
        outer_low, outer_high = max(v.lower() for v in targets), min(v.upper() for v in targets)
        inner_low, inner_high = max(a[k].upper() for a in high), min(b[k].lower() for b in low)
        def height(time, g3):
            return (coeff*(-2*time).exp()/g3-q(offset)).sqrt()
        if outer_low <= outer_high:
            possible.append((height(stop, outer_high).lower(),height(start, outer_low).upper()))
        if inner_low < inner_high:
            assert q(10**10).log() < start and stop < q(10**14).log()
            sufficient.append((height(stop, inner_high).upper(),height(start, inner_low).lower()))
    # The moving sufficient intervals cover their union as L ranges continuously over one slab.
    assert all(lo < hi for lo,hi in sufficient)
    sufficient.sort(key=lambda v: v[0])
    assert all(a[1] > b[0] for a,b in zip(sufficient,sufficient[1:]))
    outer = (min(lo for lo,hi in possible), max(hi for lo,hi in possible))
    inner = (sufficient[0][0], max(hi for lo,hi in sufficient))
    for label, ends in [('outer',outer),('inner',inner)]:
        counts = [x.zeta_nzeros().unique_fmpz() for x in ends]
        assert all(v is not None for v in counts)
        labels = (int(counts[0])+1, int(counts[1]))
        assert labels == ((10680, 12183) if label == 'outer' else (11077, 11966))
        print(label, [v.str(11) for v in ends], labels, int(counts[1]-counts[0]))


if __name__ == '__main__':
    verify()
