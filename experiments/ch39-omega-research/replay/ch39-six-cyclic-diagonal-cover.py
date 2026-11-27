"""六壳共用密度矩阵的必要条件。闭盒端点向外取整，未排除的盒始终保留。"""

from functools import lru_cache
from math import isqrt

from flint import acb_mat

N = 1 << 18
A = 1 << 32


def floor(x):
    return int(x.lower().floor().fmpz())


def ceil(x):
    return int(x.upper().ceil().fmpz())


def up(a, b):
    return -(-a // b)


def root(x, p):
    if p == 1:
        return x
    s = isqrt(x)
    return s if p == 2 else isqrt(s)


@lru_cache(maxsize=131072)
def project(h, c):
    if h * N >= c * c:
        return N
    s = isqrt(N * N - c * c) + 1
    a = isqrt(h * N) + 1
    b = isqrt((N - h) * N) + 1
    return min(N, up((c * a + s * b) ** 2, N**3))


def prepare(model):
    wl = []
    wh = []
    rows = []
    for s in model:
        for i in s["big"]:
            w = -s["w"][i]
            wl.append(floor(w * A))
            wh.append(ceil(w * A))
            rows.append([s["main"][i, k] for k in range(6)])
    peers = []
    for a, r in enumerate(rows):
        for b in range(a):
            if a // 6 == b // 6:
                continue
            c = floor(abs(sum(x * y.conjugate() for x, y in zip(r, rows[b]))) * N)
            if c > N * 6 // 10:
                peers.append((a, b, min(c, N)))
    coef = (
        -416,
        63,
        96,
        64,
        91,
        103,
        51,
        50,
        -33,
        -14,
        -24,
        -30,
        363,
        -57,
        -79,
        -37,
        -92,
        -98,
        279,
        -65,
        -51,
        -15,
        -63,
        -84,
        -1130,
        180,
        238,
        173,
        259,
        281,
        1024,
        0,
        0,
        0,
        0,
        0,
    )
    mat = acb_mat(
        [
            [
                sum(c * r[a].conjugate() * r[b] for c, r in zip(coef, rows))
                for b in range(6)
            ]
            for a in range(6)
        ]
    )
    for bound, sign in ((157, 1), (186, -1)):
        for k in range(1, 7):
            minor = acb_mat(
                [
                    [sign * (mat[i, j] - bound * int(i == j)) for j in range(k)]
                    for i in range(k)
                ]
            )
            assert minor.det().real > 0
    return wl, wh, peers, coef


def test(modeldata, targets, limit=300):
    wl, wh, peers, coef = modeldata
    conditions = []
    for j, (u, v) in enumerate(targets):
        ids = range(6 * j, 6 * j + 6)
        for p, (low, high) in (
            (1, (A, A)),
            (2, (floor(A * u), ceil(A * u))),
            (4, (floor(A * v), ceil(A * v))),
        ):
            a = [A] * 6 if p == 1 else wl[ids.start : ids.stop]
            b = [A] * 6 if p == 1 else wh[ids.start : ids.stop]
            conditions.append((ids, p, low * N**p, high * N**p, a, b))

    def tighten(l, h):
        for _ in range(12):
            prev = sum(h) - sum(l)
            for indices, p, L, U, a, b in conditions:
                fl = [aa * l[i] ** p for i, aa in zip(indices, a)]
                fh = [bb * h[i] ** p for i, bb in zip(indices, b)]
                low = sum(fl)
                high = sum(fh)
                if low > U or high < L:
                    return False
                for i, aa, bb, lc, hc in zip(indices, a, b, fl, fh):
                    l[i] = max(l[i], root(max(0, (L - high + hc) // bb), p))
                    remaining = U - low + lc
                    if remaining < 0:
                        return False
                    h[i] = min(h[i], root(up(remaining, aa), p) + 1)
                    if h[i] < l[i]:
                        return False
            for a, b, c in peers:
                l[a] = max(l[a], N - project(N - l[b], c))
                h[a] = min(h[a], project(h[b], c))
                if l[a] > h[a]:
                    return False
                l[b] = max(l[b], N - project(N - l[a], c))
                h[b] = min(h[b], project(h[a], c))
                if l[b] > h[b]:
                    return False
            # 同一迹一正算子的三十六个对角读数满足一个联合线性界。
            small = sum(c * (l[i] if c > 0 else h[i]) for i, c in enumerate(coef))
            large = sum(c * (h[i] if c > 0 else l[i]) for i, c in enumerate(coef))
            if small > 186 * N or large < 157 * N:
                return False
            for i, c in enumerate(coef):
                if c > 0:
                    lo = h[i] + (157 * N - large) // c
                    hi = l[i] + up(186 * N - small, c)
                elif c < 0:
                    lo = h[i] + (186 * N - small) // c
                    hi = l[i] + up(157 * N - large, c)
                else:
                    continue
                l[i] = max(l[i], lo)
                h[i] = min(h[i], hi)
                if l[i] > h[i]:
                    return False
            if prev - sum(h) + sum(l) < 20:
                break
        return True

    boxes = [([0] * 36, [N] * 36)]
    visits = 0
    while boxes and visits < limit:
        l, h = boxes.pop()
        visits += 1
        if not tighten(l, h):
            continue
        i = max(
            range(36),
            key=lambda k: (h[k] - l[k]) * (1 + min(10, isqrt(4 * wh[k] // A))),
        )
        if h[i] - l[i] < N // 200:
            return True, visits
        m = (h[i] + l[i]) // 2
        r = l.copy()
        r[i] = m
        u = h.copy()
        u[i] = m
        boxes.extend(((r, h), (l, u)))
    return bool(boxes), visits
