"""荷六第一算术枢轴的相邻壳投影输运，按冻结路径作单步、两步留出。"""
from pathlib import Path
from runpy import run_path

from flint import arb, ctx


CENTER = 2_306_888_100_085


def readout(number, amplitudes, weights, field, moments, heat, low_weak, high_weak):
    mp = arb('1.22089e19')
    floor = arb('91.1876')/mp
    p0, p1 = floor**2, (arb('1.176397e14')/mp)**2
    assert len(amplitudes) == len(weights) == 19
    assert all(floor < f < arb('.999') for f in amplitudes)
    assert sum(weights).contains(1)
    columns = []
    for f in amplitudes:
        x = f**2
        weak = low_weak(x/p0) if x < 1024*p0 else high_weak(x, p0, p1)
        columns.append((*moments(field, heat, x), weak))
    totals = [sum(w*v for w, v in zip(weights, col, strict=True)) for col in zip(*columns, strict=True)]
    cutoffs = [min(f.lower() for f in amplitudes), max(f.upper() for f in amplitudes)]
    goals = (arb(1), arb('1e-122'), arb(0))
    separated = [not value.overlaps(goal) for value, goal in zip(totals, goals, strict=True)]
    print('编号、截止振幅两端、K、W、D及与目标分离', number,
          *[v.str(15) for v in (*cutoffs, *totals)], separated, flush=True)


def verify(root):
    train, _ = run_path(str(root/'ch39-three-shell-two-space-readout.py'))['verify'](root, return_data=True)
    assert len(train) == 171
    lib = run_path(str(root/'ch39-adjacent-pvm-transport.py'))
    shell = run_path(str(root/'ch39-three-ray-chiral-boundary.py'))
    pivots = run_path(str(root/'ch39-native-rh-ray-pvm.py'))
    gravity = run_path(str(root/'ch39-three-index-free-atoms.py'))
    high_weak = run_path(str(root/'ch39-full-neutral-free-cutoff.py'))['weak_column']
    outer = [CENTER+k for k in (-2, -3, 2, 3)]
    moments, heat, fields = gravity['broken_fields'](root, outer)
    by_number = dict(zip(outer, fields, strict=True))
    ctx.prec = 320
    weak = pivots['weak_series']()
    for direction in (-1, 1):
        levels, previous = lib['spectrum'](shell, 6, CENTER+direction)
        amplitudes = lib['training_values'](levels, train)
        for steps in (2, 3):
            number = CENTER + direction*steps
            levels, current = lib['spectrum'](shell, 6, number)
            amplitudes = lib['transport'](previous, current, amplitudes)
            previous = current
            _, ordered, cross, _ = next(v for v in pivots['shell_weights'](root, number) if v[0] == 6)
            ctx.prec = 320
            assert all(a.overlaps(b) for a, b in zip(levels, ordered, strict=True))
            readout(number, amplitudes, [w.real for w in cross[0]], by_number[number], moments,
                    heat, weak, high_weak)


if __name__ == '__main__':
    verify(Path(__file__).resolve().parent)
