"""十九点磁壳的单步边权、两步相位重构和实际零点的有限窗口分辨率。"""
from pathlib import Path
from runpy import run_path

from flint import acb, acb_mat, arb, ctx


ROOT = Path(__file__).parent
SHELL = run_path(str(ROOT / "ch39-arc-relative-determinant.py"))
EDGES = SHELL["edges"]
CLOCKS = (724, 804)
PORTS = (1, 3)
POLE = 14


def magnetic_matrix(weights, phases):
    lap = acb_mat(19, 19)
    for (i, j, chord, native), w in zip(EDGES, weights, strict=True):
        lap[i,i] += w
        lap[j,j] += w
        lap[i,j] = -w * native * phases.get(chord, 1)
        lap[j,i] = lap[i,j].conjugate()
    return lap


def path_coefficients(weights):
    # 极点入射的全部原生相位均为一，两目标的直接边也没有零点钟。
    base = magnetic_matrix(weights, {})
    lookup = {(i,j): d for i,j,d,_ in EDGES}
    lookup.update({(j,i): d for i,j,d,_ in EDGES})
    assert all(base[k,POLE].imag.is_zero() for k in range(19))
    coeffs = acb_mat(2,2)
    offsets = []
    for row, a in enumerate(PORTS):
        constant = acb(0)
        counts = [0,0]
        for k in range(19):
            term = base[a,k]*base[k,POLE]
            chord = lookup.get((a,k))
            if chord in CLOCKS:
                assert a < k and lookup[k,POLE] not in CLOCKS
                col = CLOCKS.index(chord)
                coeffs[row,col] += term
                counts[col] += 1
            else:
                constant += term
        assert counts == [1,1]
        offsets.append(constant)
    return coeffs, acb_mat([[v] for v in offsets])


def reconstruct(truth, phases):
    lap = magnetic_matrix(truth, phases)
    second = lap*lap
    weights = [abs(lap[i,j]) for i,j,_,_ in EDGES]
    assert all(w > 0 and w.contains(v) for w,v in zip(weights,truth))
    coeffs, offsets = path_coefficients(weights)
    readings = []
    for a in PORTS:
        x, y = lap[a,POLE], second[a,POLE]
        assert x.real < 0 and x.imag.is_zero()
        # 三组强度加上已取得的 |x|²，可用极化恒等式读出 conj(x)y。
        base_intensity = abs(x)**2 + abs(y)**2
        real = (abs(x+y)**2-base_intensity)/2
        imag = (base_intensity-abs(x+acb(0,1)*y)**2)/2
        readings.append(acb(real,imag)/x.real)
    found = coeffs.solve(acb_mat([[v] for v in readings])-offsets)
    assert all(found[i,0].contains(phases[d]) for i,d in enumerate(CLOCKS))
    assert all(abs(found[i,0]-phases[d]) < arb("1e-45") for i,d in enumerate(CLOCKS))
    return coeffs


def finite_window(coeffs):
    start, count = 10855, 1238
    zeros = acb.zeta_zeros(start,count)
    observed = []
    for zero in zeros:
        z = SHELL["phases"](zero.imag)
        output = [sum(coeffs[i,j]*z[d] for j,d in enumerate(CLOCKS)) for i in range(2)]
        observed.append([part for v in output for part in (v.real,v.imag)])
    order = sorted(range(count),key=lambda i: float(observed[i][0]))
    assert all(observed[i][0] < observed[j][0] for i,j in zip(order,order[1:]))

    def distance(i,j):
        # 区间的 max，即使坐标区间重叠也给出可靠包围。
        value = arb(0)
        for x,y in zip(observed[i],observed[j]):
            value = value.max(abs(x-y))
        return value

    # 先用实部相邻的一对给出上界，再遍历此上界内的所有候选。
    upper = min((distance(i,j).upper() for i,j in zip(order,order[1:])),key=float)
    candidates = []
    for pos,i in enumerate(order):
        for j in order[pos+1:]:
            if observed[j][0]-observed[i][0] > upper:
                break
            value = distance(i,j)
            if not value > upper:
                candidates.append((value, *sorted((i+start,j+start))))
    best = min(candidates,key=lambda item:float(item[0]))
    assert best[0]>0 and all(best[0]<v for v,i,j in candidates if (i,j)!=best[1:])
    target = 11431-start
    others = [(distance(target,i),i+start) for i in range(count) if i!=target]
    nearest = min(others,key=lambda item:float(item[0]))
    assert all(nearest[0]<v for v,i in others if i!=nearest[1])
    print("全窗四个实分量的最小最大范数间隔", best[0].str(15), best[1:])
    print("第 11431 枚的最近编号和最大范数间隔", nearest[1], nearest[0].str(15))
    return best, nearest


def verify():
    ctx.prec = 192
    z = SHELL["phases"](acb.zeta_zero(11431).imag)
    one = [arb(1)]*len(EDGES)
    coeffs = reconstruct(one,z)
    varied = [arb(2+(i*31+j*17)%43)/17 for i,j,_,_ in EDGES]
    reconstruct(varied,z)
    s, t = coeffs[0,0]*coeffs[1,1], coeffs[0,1]*coeffs[1,0]
    sine_half = abs(s/abs(s)-t/abs(t))/2
    assert arb(".107") < sine_half < arb(".108")
    print("任意正边权的相对行列式下界系数", sine_half.str(15))
    best, nearest = finite_window(coeffs)
    assert best[1:] == (11282,11955) and arb(".00345810") < best[0] < arb(".00345811")
    assert nearest[1] == 11032 and arb(".1413899") < nearest[0] < arb(".1413900")

    # 与无磁第九阶读数相同的端口和内部边。L[d] 在该内部边的输运差不含钟边。
    a,b,i,j = 17,18,11,16
    support = {(min(x,y), max(x,y)):chord for x,y,chord,_ in EDGES}
    assert all(support.get(tuple(sorted((x,y)))) not in CLOCKS
               for x in (i,j) for y in (a,b))
    assert support[a,b]==56 and support[i,j]==224
    lap = magnetic_matrix(one,z)
    d = acb_mat(19,1)
    v = acb_mat(19,1)
    d[a,0], d[b,0] = 1, lap[b,a]
    v[i,0], v[j,0] = 1, lap[j,i]
    amplitude = (v.conjugate().transpose()*lap*d)[0,0]
    jump = abs(amplitude)**2/2
    assert arb("1.99") < jump < 2
    print("同一内部边加倍后第三热矩增加",jump.str(16))


if __name__ == "__main__":
    verify()
