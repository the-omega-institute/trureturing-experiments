"""归一化简单厄米铅笔本征线的区间导数与物理矩阵元。"""
from flint import arb, acb, acb_mat

def adj(x):
    return x.conjugate().transpose()

def inn(x, y):
    return (adj(x) * y)[0, 0]

def norm(x):
    return inn(x, x).real.sqrt()

def xc(x, y):
    return acb_mat([[(x[(i + 1) % 3, 0] * y[(i + 2) % 3, 0] - x[(i + 2) % 3, 0] * y[(i + 1) % 3, 0]).conjugate()] for i in range(3)])

def mats(q, remote, V):
    if remote == 1:
        return [V * q, q]
    return [q, adj(V) * q]

def centered(rootlib, geom, V, mass, A, B, chart, idx, lo, hi):
    """在中点本征矢的正相位规范中包围一阶导数，用中值定理收紧矩阵元。"""
    lib = rootlib
    H, J = (A, B) if chart == 0 else (B, A)
    remote = 1 - chart
    mid = (lo + hi) / 2
    rad = (hi - lo) / 2
    portion = lo.union(hi)
    dx = arb(0, rad.upper())
    vals, Q = (H + J * mid).eig(right=True, algorithm='rump')
    order = sorted(range(3), key=lambda i: vals[i].real.mid())
    v0 = geom['unit'](acb_mat([[Q[i, order[idx]]] for i in range(3)]))
    ls = lib['lines'](H, J, lo, hi, (0, 1, 2), geom)
    if len(ls) < 3:
        return None
    eig = [vals[i].real + dx for i in order]
    if any((eig[idx].overlaps(eig[j]) for j in range(3) if j != idx)):
        return None
    y = ls[idx]
    p = inn(v0, y)
    if not abs(p) > 0:
        return None
    y = y * (p.conjugate() / abs(p))
    jp = J * y - y * inn(y, J * y).real.mid()
    dy = acb_mat(3, 1)
    for i in range(3):
        if i != idx:
            dy += ls[i] * (inn(ls[i], jp) / (eig[idx] - eig[i]))
    overlap = inn(v0, y).real
    if not overlap > 0:
        return None
    dy += y * acb(0, -inn(v0, dy).imag / overlap)
    diag = [geom['diag'](v) for v in mass]

    def frame(yy, deriv=None):
        yrep = V * yy if remote == 0 else yy
        qr = diag[remote]
        j = max(range(3), key=lambda j: abs(yrep[j, 0]).mid())
        weight = mass[remote][j]
        sy = weight + sum(((mass[remote][k] - weight) * abs(yrep[k, 0]) ** 2 for k in range(3)))
        g = acb_mat([[(mass[remote][i] - sy) * yrep[i, 0]] for i in range(3)])
        w = norm(g)
        if not w > 0:
            return None
        zz = g * (1 / w)
        xx = xc(yrep, zz)
        reps = [mats(q, remote, V) for q in (xx, yrep, zz)]
        if deriv is not None:
            yp = V * deriv if remote == 0 else deriv
            sp = 2 * sum(((mass[remote][k] - weight) * (yp[k, 0].conjugate() * yrep[k, 0]).real for k in range(3)))
            gp = acb_mat([[(mass[remote][i] - sy) * yp[i, 0] - yrep[i, 0] * sp] for i in range(3)])
            wp = inn(zz, gp).real
            zp = (gp - zz * wp) * (1 / w)
            xp = xc(yp, zz) + xc(yrep, zp)
            derivs = [mats(q, remote, V) for q in (xp, yp, zp)]
        out = []
        for i, D in enumerate(diag):
            xh, yh, zh = [r[i] for r in reps]

            def shifted(h):
                jj = max(range(3), key=lambda j: abs(h[j, 0]).mid())
                mm = mass[i][jj]
                return (mm, acb_mat([[(mass[i][k] - mm) * h[k, 0]] for k in range(3)]))
            mx, X = shifted(xh)
            my, Y = shifted(yh)
            mz, Z = shifted(zh)
            K = D * zh - zh * mx
            W = acb_mat([[(mass[i][k] - my) * zh[k, 0]] for k in range(3)])
            if deriv is None:
                out.append((mx + inn(xh, X).real, my + inn(yh, Y).real, mz + inn(zh, Z).real, inn(xh, K), inn(yh, W).real))
            else:
                xp, yp, zp = [r[i] for r in derivs]
                out.append((2 * inn(xp, X).real, 2 * inn(yp, Y).real, 2 * inn(zp, Z).real, inn(xp, K) + inn(X, zp), (inn(yp, W) + inn(Y, zp)).real))
        return out
    val = frame(v0)
    deriv = frame(y, dy)
    if val is None or deriv is None:
        return None
    out = []
    for i, (v0s, dd) in enumerate(zip(val, deriv)):
        arr = [z + dx * d for z, d in zip(v0s, dd)]
        if not arr[0] > 0 or not arr[1] > 0 or (not 1 - arr[0] > 0) or (not 1 - arr[1] > 0):
            return None
        out.append((*arr, mass[i][0] * mass[i][1] * mass[i][2]))
    return out
