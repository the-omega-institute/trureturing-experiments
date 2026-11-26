"""Certified discrepancies of deterministic zero-to-corner assignments on the shortest triangle."""
import json
from fractions import Fraction as F
from itertools import combinations
from math import ceil, floor
from pathlib import Path
from runpy import run_path

ROOT = Path(__file__).resolve().parent
THIRD = F(1, 3)


def raw_weights():
    records = json.loads((ROOT/'ch39-tick-separators.json').read_text())['prefix']
    raw = []
    unit = 1 << 100
    for cell, *_ in records[:16]:
        lo, hi = F(cell, 1 << 400), F(cell+1, 1 << 400)
        wlo, whi = 1/(F(1, 4)+hi**2), 1/(F(1, 4)+lo**2)
        raw.append((F(floor(wlo*unit), unit), F(ceil(whi*unit), unit)))
    return raw


def shares(raw, indices):
    lower = sum(raw[i][0] for i in indices)
    upper = sum(raw[i][1] for i in indices)
    return lower/sum(b for a, b in raw), upper/sum(a for a, b in raw)


def certify_fifteen(raw, second, deviation_upper):
    # Index zero stands alone. For every partition other than the displayed
    # pair, one side is already strictly below the displayed worst share.
    assert shares(raw, [0, 14])[0] > THIRD+deviation_upper
    unit = 1 << 100
    upper = [int(b*unit) for _, b in raw[1:]]
    lower_total = sum(a for a, _ in raw)*unit
    cutoff = (THIRD-deviation_upper)*lower_total
    subset = [0]
    for weight in upper:
        subset += [s+weight for s in subset]
    witness = sum(1 << (i-1) for i in second)
    full = len(subset)-1
    for mask in range(len(subset)):
        if mask not in (witness, full ^ witness):
            assert min(subset[mask], subset[full ^ mask]) < cutoff


def verify():
    tri = run_path(str(ROOT/'ch39-weyl-triangle-average.py'))
    x, z, den = tri['X'], tri['Z'], tri['denominator']
    y = tri['third_vertex'](den(x, z)/(3*den(x, x)-den(x, z)))
    contrast = [a-b for a, b in zip(tri['edge'](x, x), tri['edge'](x, y), strict=True)]
    span = max(contrast)-min(contrast)
    reference = tri['edge'](x, z)[0]
    raw = raw_weights()
    assert all(x[0] > y[1] > 0 for x, y in zip(raw[:-1], raw[1:], strict=True))
    for n, second, lower, upper in [
        (7, [2, 3], '0.07759598427', '0.07759598428'),
        (15, [2, 5, 6, 9, 10, 14], '0.00004918553', '0.00004918554'),
        (16, [2, 5, 7, 9, 11, 13, 15], '0.00483274468', '0.00483274469'),
    ]:
        prefix = raw[:n]
        lo, hi = shares(prefix, [0])
        second = {i-1 for i in second}
        third = set(range(1, n))-second
        if n == 15:
            l3, u3 = shares(prefix, third)
            a, b = THIRD-u3, THIRD-l3
            assert hi-THIRD < a
            certify_fifteen(prefix, second, b)
        elif lo > THIRD:
            a, b = lo-THIRD, hi-THIRD
        else:
            a, b = THIRD-hi, THIRD-lo
            assert lo+shares(prefix, [n-1])[0] > THIRD+b
        assert F(lower) < a < b < F(upper)
        for group in ((second,) if n == 15 else (second, third)):
            gl, gu = shares(prefix, group)
            assert THIRD-a < gl < gu < THIRD+a
        print(n, 'sharp share deviation in', lower, upper)
        if n == 7:
            assert F('0.11334308') < a*span/reference < b*span/reference < F('0.11334309')
            assert F('0.18883048') < a/hi < b/lo < F('0.18883050')
            for size in range(n+1):
                for group in combinations(range(n), size):
                    gl, gu = shares(prefix, group)
                    assert gu < THIRD-F('0.00047292') or gl > THIRD+F('0.00047292')
            # Fill the last two bins by splitting the first atom across all three.
            assert shares(prefix, second)[1] < THIRD
            assert shares(prefix, set(range(1, n))-second)[1] < THIRD
    assert shares(raw[:14], [0])[0] > THIRD+F(1, 1000)
    print('First prefix with a share discrepancy below 0.0001: 15. Seven atoms need nine placements.')


if __name__ == '__main__':
    verify()
