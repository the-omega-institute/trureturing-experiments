"""第四步已给出第五壳正解时，第六壳任意程度的非选择性测量不能恢复第七壳真空。"""

from pathlib import Path
from runpy import run_path

from flint import arb, arb_mat, ctx


def transition(target, source):
    amplitudes = target.conjugate().transpose()*source
    return arb_mat([[abs(a)**2 for a in row] for row in amplitudes.tolist()])


def verify(root):
    def lib(name):
        return run_path(str(root/("ch39-"+name+".py")))

    fifth, sixth, seventh = (2_306_888_100_085 + j for j in (5, 6, 7))
    spectrum = lib("adjacent-pvm-transport")["spectrum"]
    shell = lib("three-ray-chiral-boundary")
    b5, b6, b7 = (spectrum(shell, 6, n)[1] for n in (fifth, sixth, seventh))
    direct = transition(b7, b5)
    via = transition(b7, b6)*transition(b6, b5)
    assert all(sum(row).contains(1) for matrix in (direct, via) for row in matrix.tolist())
    _, ordered = lib("adjacent-ordered-charge-six-cutoffs")["verify"](
        root, refine_node=2_306_888_100_086)[2]
    weights = []
    for node in (fifth, seventh):
        _, _, pivots, _ = next(row for row in lib("native-rh-ray-pvm")["shell_weights"](root,node)
                               if row[0] == 6)
        weights.append([v.real for v in pivots[2]])
    _, heat, fields = lib("three-index-free-atoms")["broken_fields"](root, [fifth, seventh])
    ctx.prec = 192
    c, e = weights
    trace = sum(v.sqrt() for v in ordered)
    strong = {0, 9, 15, 16, 17, 18}
    low = [arb(0) for _ in c]
    low[2], low[3], low[14] = (arb(v).sqrt().lower() for v in (".025",".10",".001"))
    caps = [".001" if j in strong else ".75553" for j in range(19)]
    assert ordered[-1].sqrt() < arb(".75553")
    assert [i+1 for i, value in enumerate(e) if value>0] == [1, 7, 10, 16]
    # 任一中间退相干程度给出逐行随机矩阵。正权行作凸上界，负权行作下界。
    delta = arb(".04")
    coeff_positive = [arb(0) for _ in range(19)]
    negative = [[] for _ in range(19)]
    shift = arb(0)
    for i, weight in enumerate(e):
        row = [direct[i,j].union(via[i,j]) for j in range(19)]
        if weight > 0 and i in strong:
            cap = sum(value*arb(limit) for value,limit in zip(row,caps,strict=True))
            shift += weight*(cap.upper()+delta)**4
        elif weight > 0:
            for j, value in enumerate(row):
                coeff_positive[j] += weight*value.upper()
        else:
            assert weight < 0
            for j, value in enumerate(row):
                negative[j].append((weight, max(arb(0),value.lower())))
    # 七壳振幅可与预测读数逐层相差 delta。对负权行利用
    # Σ(p_j y_j−delta)_+ ≤ (Σp_j y_j−delta)_+，再用四次幂的超可加性。
    u, v, d = map(arb, (".47","-1.98",".162"))
    errors = lib("fifth-seventh-cutoff-separation")["moment_errors"]
    s2, s4 = errors(fields[0],heat,c,caps)
    base = shift + u*s2 + v*s4 + d*trace
    parts = 4096
    for i in range(19):
        lo, hi = low[i], arb(caps[i])
        peak = None
        for step in range(parts):
            left = lo + (hi-lo)*step/parts
            right = lo + (hi-lo)*(step+1)/parts
            square = (left**2).union(right**2)
            high = coeff_positive[i]*(right+delta)**4
            high += sum(weight*max(arb(0),(p*left-delta).lower())**4
                        for weight,p in negative[i])
            bound = (high-v*c[i]*square**2-u*c[i]*square-d*left).upper()
            peak = bound if peak is None else max(peak,bound)
        base += peak
    # 真空热矩 W=w*S4+R。用无质量主项和包围余项将 S4 的界还原为 W。
    _, target = errors(fields[1],heat,e,[".75553"]*19)
    scale = -sum(g for g,_,_ in fields[1])/(8*arb.pi())
    vacuum = arb("1e-122") + scale*(base-target)
    assert scale>0
    print("去相干预测的 0.04 邻域内，第七壳四次矩上界、实际真空项上界",
          base.upper().str(16),vacuum.upper().str(16),flush=True)
    assert vacuum < arb("-.19769")
    return vacuum.upper()


if __name__ == "__main__":
    verify(Path(__file__).resolve().parent)
