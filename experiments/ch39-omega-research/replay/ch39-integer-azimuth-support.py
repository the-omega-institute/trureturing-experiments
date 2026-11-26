"""将原生下型支撑不等式的完整方位角覆盖化为定点整数运算。"""
from ctypes import c_int64, c_int, POINTER, CDLL
from functools import lru_cache
from pathlib import Path
from platform import system
from shutil import which
from subprocess import run
from tempfile import TemporaryDirectory

import numpy as np
from flint import arb

M = 10000000
ptr = POINTER(c_int64)
support_library = None

def compile_cover(root):
    global support_library
    memos.clear()
    cache.clear()
    cc = which('cc')
    if cc is None:
        raise RuntimeError('精确整数方位角检验需要支持 128 位整数的 C 编译器 cc')
    with TemporaryDirectory() as temp:
        output = Path(temp) / 'azimuth.so'
        option = '-dynamiclib' if system() == 'Darwin' else '-shared'
        run([cc, '-std=c99', '-O3', '-fPIC', option, str(root / 'ch39-azimuth-integer-cover.c'), '-o', str(output)], check=True)
        support_library = CDLL(str(output))
    support_library.supports.argtypes = [c_int64, c_int64, ptr, ptr, ptr]
    support_library.supports.restype = c_int

def lower(x):
    return int((M * x).lower().floor().unique_fmpz())

def upper(x):
    return int((M * x).upper().ceil().unique_fmpz())

def point(a):
    return a.ctypes.data_as(ptr)

@lru_cache(maxsize=1536)
def circle(depth, j):
    return (arb(2 * j + 1, 1) * arb.pi() / 2 ** depth).sin_cos()
memos = {}
cache = {}

def initialize(forms):
    idf = id(forms)
    if idf not in memos:
        flat = np.array([[lower(b[0]), upper(b[1]), lower(b[2]), upper(b[2]), upper(b[5])] for b in forms], dtype=np.int64)
        assert np.max(np.abs(flat)) < 8 * M
        memos[idf] = (flat, {})
    return memos[idf]

def section(forms, depth, ranges):
    if depth not in ranges:
        vals = []
        for j in range(2 ** depth):
            sn, cs = circle(depth, j)
            l = []
            for n0, sup, a, b, c, n2 in forms:
                v = b * cs + c * sn
                l.append((lower(v), upper(v)))
            vals.append(l)
        ranges[depth] = np.array(vals, dtype=np.int64)
        assert np.max(np.abs(ranges[depth])) < 4 * M
    return ranges[depth]

def test(forms, box):
    """将每段的三角区间与支撑系数向外取整，再检查全部方位角段。"""
    if any(abs(value) > 4 * M for value in box):
        return False
    flat, ranges = initialize(forms)
    key = (id(forms), *box)
    if key in cache:
        return cache[key]
    dat = np.array(box, dtype=np.int64)
    good = False
    for depth in (4, 6, 8):
        interval = section(forms, depth, ranges)
        if support_library.supports(len(forms), 2 ** depth, point(flat), point(interval), point(dat)):
            good = True
            break
    if len(cache) < 1000000:
        cache[key] = good
    return good

def upper_grid(n, m):
    return -(-n // m) * m

def bound_orientation(row, inv):
    j, cap, dist, forms = row
    tu, td, r, p, l, polar, rd = inv
    de = 1 - r ** 2
    if dist * (1 - cap / tu) > de:
        return True
    if not r > 0:
        return False
    K = 2 / dist
    if not K * de < 1:
        return False
    h = (1 - K * de).sqrt()
    E = ((K - 1) * de ** 2 / (r * (r + h) * (1 + h) * (1 + r))).upper()
    si = (2 * E).sqrt()
    if not all((v.is_finite() for v in (p, l, td, rd, E, si))):
        return False
    box = (lower(p), upper(p), lower(l), upper(l), lower(td), upper(rd), upper(E), upper(si))
    pp, ph, ll, lh, tl, rhi, ef, si = box
    for grid in (30000, 2500):
        grain = max(1, tl // (400 if grid == 30000 else 4000))
        q = (pp // grid * grid, upper_grid(ph, grid), ll // grid * grid, upper_grid(lh, grid), tl // grain * grain, upper_grid(rhi, grid), upper_grid(ef, grid), upper_grid(si, grid))
        if test(forms, q):
            return True
    return test(forms, box)
